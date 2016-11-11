import unittest
import random
from driven.generic.adapter import full_genotype, ModelModificationMixin, GenotypeChangeModel, \
    MediumChangeModel, MeasurementChangeModel
from cameo import models
from cobra.manipulation.delete import find_gene_knockout_reactions


GENES_TO_ADD = [
    'DDC',
    'tpH',
    'marC',
    'TpH_mN',
    'asmt',
    'AANAT',
    'FolE',
    'PhhB',
    'CytR_delta8bp',
    'ASMT',
    'TsgA',
    'PtrpE',
    'PJ23119',
    'TrpE',
    'kanMx',
    'ptrc',
    'argG_mA',
    'Cat-SacB',
    'sgAANAT',
    'pphB',
]

REACTIONS_TO_EQUATIONS = {
    'ko:K01724': 'C00002 + C01300 <=> C00020 + C04807',
    'ko:K00543': 'C00019 + C05635 <=> C00021 + C05660',
    'ko:K00666': 'C00019 + C00002 <=> C00020 + C01300',
}

GENES_TO_REACTIONS = {
    key: REACTIONS_TO_EQUATIONS for key in GENES_TO_ADD
}


# ModelModificationMixin is a mixin, instances shoudn't be created
class TestAdapter(ModelModificationMixin):
    def __init__(self, model):
        self.added_reactions = set()
        self.model = model


class AdapterTestCase(unittest.TestCase):
    def test_model_modification(self):
        model = models.bigg.iJO1366
        adapter = TestAdapter(model)
        metabolite = model.metabolites.get_by_id('13dpg_c')
        self.assertNotIn('13dpg_e', model.metabolites)
        adapter.create_exchange(metabolite)
        self.assertIn('13dpg_e', model.metabolites)
        self.assertEqual(model.metabolites.get_by_id('13dpg_e').formula, metabolite.formula)


def random_genes_to_knockout(model):
    genes_to_knock = random.sample(model.genes, 10)
    return {gene.name: gene.id for gene in genes_to_knock}


def check_knockout_bounds(model, gene, function):
    for reaction in find_gene_knockout_reactions(model, [gene]):
        function(reaction.lower_bound == reaction.upper_bound == 0)

wild_model = models.bigg.iJO1366


class TestStrainToModel(unittest.TestCase):
    def test_knockout_reaction_bounds(self):
        """If gene is knocked out correctly, reaction's lower and upper bounds should be set to zero"""
        model = wild_model.copy()
        check_knockout_bounds(model, model.genes.b2103, self.assertFalse)
        model = GenotypeChangeModel(wild_model.copy(), full_genotype(['-thiD']), GENES_TO_REACTIONS).model
        check_knockout_bounds(model, model.genes.b2103, self.assertTrue)

    def test_chains_of_knockouts(self):
        """Chains of knockouts should be performed correctly"""
        model = wild_model.copy()
        genes = random_genes_to_knockout(model)
        knockout_chain = ['-' + g for g in genes.keys()]
        model = GenotypeChangeModel(wild_model.copy(), full_genotype(knockout_chain), GENES_TO_REACTIONS)
        for gene_name, gene_id in genes.items():
            check_knockout_bounds(model.model, getattr(model.model.genes, gene_id), self.assertTrue)
        self.assertTrue(set(genes.keys()) == model.knocked_out_genes)

    def test_add_genes(self):
        """Adding a gene makes the new reactions appear"""
        gene_name = 'pphB'
        for reaction in GENES_TO_REACTIONS[gene_name]:
            model = GenotypeChangeModel(wild_model.copy(), full_genotype(['+' + gene_name]), GENES_TO_REACTIONS)
            self.assertTrue(hasattr(model.model.reactions, reaction))

    def test_measurements(self):
        """Check if bounds for exchange reactions are correct"""
        model = wild_model.copy()
        measurements = [{'id': 'bigg:thrp', 'measurement': 15}, {'id': 'chebi:16695', 'measurement': -3}]
        model = MeasurementChangeModel(model, measurements).model
        self.assertEqual(model.reactions.EX_thrp_e.lower_bound, 15)
        self.assertEqual(model.reactions.EX_thrp_e.upper_bound, 15)
        self.assertEqual(model.reactions.EX_ump_e.lower_bound, -3)
        self.assertEqual(model.reactions.EX_ump_e.upper_bound, -3)

    def test_medium(self):
        model = wild_model.copy()
        """Demand reaction bounds are changed depending on carbon in the formula"""
        medium = [{'id': 'bigg:23dappa', 'concentration': 0.6}, {'id': 'chebi:23335', 'concentration': 0.05}]
        model = MediumChangeModel(model, medium).model
        self.assertEqual(model.reactions.EX_23dappa_e.lower_bound, -1)
        self.assertEqual(model.reactions.EX_23dappa_e.upper_bound, wild_model.reactions.EX_23dappa_e.upper_bound)
        self.assertEqual(model.reactions.EX_cobalt2_e.lower_bound, -1000)
        self.assertEqual(model.reactions.EX_cobalt2_e.upper_bound, wild_model.reactions.EX_cobalt2_e.upper_bound)
