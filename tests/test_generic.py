import unittest
from driven.generic.adapter import ModelModificationMixin
from cameo import models


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
