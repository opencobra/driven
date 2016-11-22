import re
import gnomic
from cameo.data import metanetx
from cameo.core.reaction import Reaction
from cameo.core.metabolite import Metabolite
import logging
logger = logging.getLogger(__name__)


def clean_bigg_id(string):
    return re.sub(r"bigg:|dsh", "", string)


def get_existing_metabolite(mnx_id, model, compartment):
    """Find compartment in the model by Metanetx id.

    Parameters
    ----------
    mnx_id : string
        Metanetx id
    model
        cobra model
    compartment : string
        f.e "_c"

    Returns
    -------
    Metabolite or None

    """
    if not mnx_id:
        return
    try:
        clean_id = clean_bigg_id(metanetx.mnx2bigg[mnx_id])
        return model.metabolites.get_by_id(clean_id + compartment)
    except KeyError:
        try:
            return model.metabolites.get_by_id(mnx_id + compartment)
        except KeyError:
            pass


def contains_carbon(metabolite):  # TODO: use method from Metabolite class when this change is merged
    if not metabolite.formula:
        raise ValueError("No formula for metabolite {}, it's unknown if there is carbon in it")
    return 'C' in metabolite.elements


def find_metabolite_info(met_id):
    """Find chemical formula of metabolite in metanetx.chem_prop dictionary

    Parameters
    ----------
    met_id : string
        string of format "<metabolite_id>_<compartment_id>", where <metabolite_id> is a BIGG id or a Metanetx id

    Returns
    -------
    pandas row or None

    """
    met_id = met_id[:-2]
    try:
        if met_id in metanetx.chem_prop.index:
            return metanetx.chem_prop.loc[met_id]
        return metanetx.chem_prop.loc[metanetx.all2mnx['bigg:' + met_id]]
    except KeyError:
        return None


def feature_id(feature):
    return feature.name if feature.name else feature.accession.identifier


class ModelModificationMixin(object):
    """
    Base model modification class, providing methods for adding adapter and exchange reactions for new metabolites
    """
    model = None
    changes = None

    def create_exchange(self, metabolite):
        """For given metabolite A_c from c compartment, create:
        a) corresponding metabolite A_e from e compartment;
        b) adapter reaction A_c <--> A_e
        c) exchange reaction A_e -->

        Parameters
        ----------
        metabolite : Metabolite
            metabolite id in format <bigg_id>_c, f.e. Nacsertn_c

        Returns
        -------

        """
        exchange_metabolite = Metabolite(metabolite.id.replace('_c', '_e'), formula=metabolite.formula, compartment='e')
        self.add_adapter_reaction(metabolite, exchange_metabolite)
        self.add_exchange_reaction(exchange_metabolite)

    def add_exchange_reaction(self, metabolite):
        """Add exchange reaction "A --> " for given metabolite A
         If reaction exists, log and pass

        Parameters
        ----------
        metabolite : basestring
            metabolite id in format <bigg_id>_<compartment_id>, f.e. Nacsertn_c

        Returns
        -------

        """
        try:
            logger.debug('Add exchange reaction for metabolite: {}'.format(metabolite.id))
            exchange_reaction = self.model.add_exchange(metabolite, prefix='EX_')
            self.changes['added']['reactions'].add(exchange_reaction)
            self.annotate_new_metabolites(exchange_reaction)
        except ValueError:
            logger.debug('Exchange reaction exists for metabolite {}'.format(metabolite.id))

    def add_adapter_reaction(self, metabolite, existing_metabolite):
        """Add adapter reaction A <--> B for metabolites A and B

        Parameters
        ----------
        metabolite : Metabolite
            metabolite A
        existing_metabolite : Metabolite
            metabolite B

        Returns
        -------

        """
        try:
            adapter_reaction = Reaction(str('adapter_' + metabolite.id + '_' + existing_metabolite.id))
            adapter_reaction.lower_bound = -1000
            adapter_reaction.add_metabolites({metabolite: -1, existing_metabolite: 1})
            self.model.add_reactions([adapter_reaction])
            self.changes['added']['reactions'].add(adapter_reaction)
            self.annotate_new_metabolites(adapter_reaction)
            logger.debug('Adapter reaction added: {} <--> {}'.format(metabolite.id, existing_metabolite.id))
        except Exception:  # TODO: raise a reasonable exception on cameo side if the reaction exists
            logger.debug('Adapter reaction exists: {} <--> {}'.format(metabolite.id, existing_metabolite.id))

    def add_demand_reaction(self, metabolite):
        """For metabolite in e compartment with existing exchange reaction, make it possible to consume metabolite
        by decreasing the lower bound of exchange reaction

        Parameters
        ----------
        metabolite : Metabolite
            metabolite from e compartment, f.e. melatn_e

        Returns
        -------

        """
        exchange_reaction = list(set(metabolite.reactions).intersection(self.model.exchanges))[0]
        if exchange_reaction.lower_bound >= 0:
            exchange_reaction.lower_bound = -1 if contains_carbon(metabolite) else -1000
        self.changes['added']['reactions'].add(exchange_reaction)
        self.annotate_new_metabolites(exchange_reaction)

    @staticmethod
    def annotate_new_metabolite(metabolite):
        """Find information about new metabolite in chem_prop dictionary and add it to model

        Parameters
        ----------
        metabolite : Metabolite
            new metabolite

        Returns
        -------

        """
        info = find_metabolite_info(metabolite.id)
        if info is not None:
            metabolite.formula = info['formula']
            metabolite.name = info['name']
            metabolite.annotation = info.to_dict()
        else:
            logger.debug('No formula for {}'.format(metabolite.id))

    def annotate_new_metabolites(self, reaction):
        """Annotate new metabolites with chem_prop information and keep track of them

        Parameters
        ----------
        reaction : Reaction
            reaction that is added to the model

        Returns
        -------

        """
        for metabolite in reaction.metabolites:
            if metabolite.formula is None:  # unknown metabolite
                self.annotate_new_metabolite(metabolite)
                self.changes['added']['metabolites'].add(metabolite)

    def model_metabolite(self, metabolite_id, compartment='_e'):
        """Get metabolite associated with this model for a given entity

        Parameters
        ----------
        metabolite_id
            string of format <database>:<id>, f.e. chebi:12345
        compartment
            the compartment where to find the metabolite, e.g. _e for exchange
        Returns
        -------
        the model metabolite (or None if no matching found)
        """
        mnx_id = metanetx.all2mnx.get(metabolite_id)
        return get_existing_metabolite(mnx_id, self.model, compartment)


def map_equation_to_bigg(equation, compartment=None):
    """Try to map given equation which contains KEGG ids to the equation which contains BIGG ids.
    If metabolite does not exist in the BIGG database, use Metanetx id.
    If compartment is given, metabolites ids will have it as postfix.

    Example:
    Input: C00002 + C00033 <=> C00013 + C05993, compartment='_c'
    Output: atp_c + ac_c <=> ppi_c + MNXM4377_c

    :param equation: string
    :param compartment: f.e. "_c"
    :return:
    """
    array = equation.split()
    result = []
    for i, el in enumerate(array):
        if not re.match("^[A-Za-z][A-Za-z0-9]*$", el):
            result.append(el)
        else:
            try:
                el = metanetx.all2mnx['kegg:' + el]
                el = metanetx.mnx2bigg[el].replace('bigg:', '')
            except KeyError:
                pass
            if compartment:
                el += compartment
            result.append(el)
    return ' '.join(result)


def full_genotype(genotype_changes):
    """Construct gnomic Genotype object from the list of strings with changes

    :param genotype_changes: list of changes, f.e. ['-tyrA::kanMX+', 'kanMX-']
    :return:
    """

    def chain(definitions, **kwargs):
        if not definitions:
            return gnomic.Genotype([])
        genotype = gnomic.Genotype.parse(definitions[0], **kwargs)
        for definition in definitions[1:]:
            genotype = gnomic.Genotype.parse(definition, parent=genotype, **kwargs)
        return genotype

    return chain(genotype_changes)


class GenotypeChangeModel(ModelModificationMixin):
    """
    Applies genotype change on cameo model
    """

    def __init__(self, model, genotype_changes, genes_to_reactions):
        """Initialize change model

        :param model: cameo model
        :param genotype_changes: gnomic.Genotype object
        :param genes_to_reactions: dictionary like {<gene name>: {<reaction id>: <reactions equation>, ...}, ...}
        """
        self.compartment = '_c'
        self.model = model
        self.genes_to_reactions = genes_to_reactions
        self.changes = {
            'added': {'reactions': set(), 'metabolites': set()},  # reaction contain information about genes
            'removed': {'genes': set()},
        }
        self.apply_changes(genotype_changes)

    def apply_changes(self, genotype_changes):
        """Apply genotype changes on initial model

        :param genotype_changes: gnomic.Genotype
        :return:
        """
        for change in genotype_changes.changes():
            if isinstance(change, gnomic.Mutation):
                self.apply_mutation(change)
            if isinstance(change, gnomic.Plasmid):
                self.add_plasmid(change)

    def apply_mutation(self, mutation):
        """Apply mutations on initial model

        :param mutation: gnomic.Mutation
        :return:
        """
        if mutation.old:
            for feature in mutation.old.features():
                self.knockout_gene(feature)
        if mutation.new:
            for feature in mutation.new.features():
                self.add_gene(feature)

    def add_plasmid(self, plasmid):
        """Add plasmid features to the initial model.
        No plasmid instance in cameo, so changes are made in model genes and reactions directly

        :param plasmid: gnomic.Plasmid
        :return:
        """
        for feature in plasmid.features():
            self.add_gene(feature)

    def knockout_gene(self, feature):
        """Perform gene knockout.
        Use feature name as gene name

        :param feature: gnomic.Feature
        :return:
        """
        gene = self.model.genes.query(feature.name, attribute="name")
        if gene:
            gene[0].knock_out()
            self.changes['removed']['genes'].add(gene[0])
            logger.info('Gene knockout: {}'.format(gene[0].name))
        else:
            logger.info('Gene for knockout is not found: {}'.format(feature.name))

    def add_gene(self, feature):
        """Perform gene insertion.
        Find all the reactions associated with this gene using KEGGClient and add them to the model

        :param feature: gnomic.Feature
        :return:
        """
        logger.info('Add gene: {}'.format(feature.name))
        identifier = feature_id(feature)
        if self.model.genes.query(identifier, attribute='name'):  # do not add if gene is already there
            logger.info('Gene {} exists in the model'.format(feature.name))
            return
        for reaction_id, equation in self.genes_to_reactions.get(identifier, {}).items():
            self.add_reaction(reaction_id, equation, identifier)
        logger.info('Gene added: {}'.format(identifier))

    def add_reaction(self, reaction_id, equation, gene_name):
        """Add new reaction by rn ID from equation, where metabolites defined by kegg ids.

        :param reaction_id: reaction rn ID
        :param equation: equation string, where metabolites are defined by kegg ids
        :param gene_name: gene name
        :return:
        """
        reaction = Reaction(reaction_id)
        self.model.add_reactions([reaction])
        equation = map_equation_to_bigg(equation, self.compartment)
        logger.info('New reaction: {}'.format(equation))
        reaction.build_reaction_from_string(equation)
        for metabolite in reaction.metabolites:
            if metabolite.formula is None:  # unknown metabolite
                self.create_exchange(metabolite)
        reaction.gene_reaction_rule = gene_name
        self.changes['added']['reactions'].add(reaction)


class MediumChangeModel(ModelModificationMixin):
    """
    Applies medium on cameo model
    """

    def __init__(self, model, medium):
        """
        Parameters
        ----------
        model
            cameo model
        medium
            list of dictionaries of format
            {'id': <compound id (<database>:<id>, f.e. chebi:12345)>, 'concentration': <compound concentration (float)>}
        """
        self.medium = medium
        self.model = model
        self.changes = {
            'added': {'reactions': set()},
        }
        self.apply_medium()

    def apply_medium(self):
        """For each metabolite in medium try to find corresponding metabolite in e compartment of the model.
        If metabolite is found, change the lower limit of the reaction to a negative number,
        so the model would be able to consume this compound.
        If metabolite is not found in e compartment, log and continue.
        """
        for compound in self.medium:
            existing_metabolite = self.model_metabolite(compound['id'], '_e')
            if not existing_metabolite:
                logger.info('No metabolite {}'.format(compound['id']))
                continue
            logger.info('Found metabolite {}'.format(compound['id']))
            self.add_demand_reaction(existing_metabolite)


class MeasurementChangeModel(ModelModificationMixin):
    """
    Update constraints based on measured fluxes
    """

    def __init__(self, model, measurements):
        """

        Parameters
        ----------
        model
            cameo model
        measurements
            list of dictionaries of format
            {'id': <metabolite id (<database>:<id>, f.e. chebi:12345)>, 'measurement': <measurement (float)>}
        """
        self.measurements = measurements
        self.model = model
        self.changes = {
            'added': {'reactions': set()},
        }
        self.missing_in_model = []
        self.apply_exchanges()

    def apply_exchanges(self):
        """For each measured flux (production-rate / uptake-rate), constrain the model by setting
        upper and lower bound locked to these values. """
        for scalar in self.measurements:
            model_metabolite = self.model_metabolite(scalar['id'], '_e')
            if not model_metabolite:
                self.missing_in_model.append(scalar['id'])
                logger.info('Model is missing metabolite {}'.format(scalar['id']))
                return
            reaction = list(set(model_metabolite.reactions).intersection(self.model.exchanges))[0]
            reaction.change_bounds(lb=scalar['measurement'], ub=scalar['measurement'])
            self.changes['added']['reactions'].add(reaction)
