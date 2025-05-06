# Porthmeus
# 24.02.23
# coding: utf-8

import unittest
from pathlib import Path
from src.MeMoModel import MeMoModel
from src.MeMoMetabolite import MeMoMetabolite


class Test_MeMoMetabolite(unittest.TestCase):

    def test_emptyInchiStringHandling(self):
        b = MeMoModel()
        metabolite: MeMoMetabolite = MeMoMetabolite()
        metabolite.set_id("A")
        metabolite.set_inchi_string("", source="test")
        metabolite.add_inchi_string("", source="test")
        metabolite.set_names(["Peter"], source="test")
        self.assertTrue(metabolite._inchi_string == None)

    def test_AnnotateMatchingOnEmptyInchis(self):
        b = MeMoModel()

        metabolite: MeMoMetabolite = MeMoMetabolite()
        metabolite.set_id("A")
        metabolite.set_inchi_string("", source="test")
        metabolite.add_inchi_string("", source="test")
        metabolite.set_names(["Peter"], source="test")

        metaboliteB: MeMoMetabolite = MeMoMetabolite()
        metaboliteB.set_id("B")
        metaboliteB.set_inchi_string("", source="test")
        metaboliteB.add_inchi_string("", source="test")
        metaboliteB.set_names(["Peter"], source="test")

        model1 = MeMoModel([metabolite])
        model2 = MeMoModel([metaboliteB])

        model1.annotate()
        model2.annotate()

        try:
            model1.match(model2)
        except KeyError:
            self.fail("This call should not fail")

    def test_sourceAnnotation(self):

        # test automatic addition of "model" with initiation
        met = MeMoMetabolite(
            _id="A",
            _inchi_string="test",
            annotations={"test": ["test"]},
            names=["Peter"],
        )

        self.assertEqual(met._inchi_source, "model")
        self.assertEqual(met.names_source, ["model"])
        self.assertEqual(met.annotations_source, {"test": ["model"]})

        # test manual addition with initiation
        met2 = MeMoMetabolite(
            _id="A",
            _inchi_string="test",
            _inchi_source="manual",
            annotations={"test": ["manual_test"]},
            annotations_source="manual",
            names=["Hans"],
            names_source="manual",
        )

        self.assertEqual(met2._inchi_source, "manual")
        self.assertEqual(met2.names_source, ["manual"])
        self.assertEqual(met2.annotations_source, {"test": ["manual"]})

        # test merging - this tests on the one hand that the correct source are carried over, as well as that the sorting is working
        met.merge(met2)
        self.assertEqual(met.names_source, ["manual", "model"])
        self.assertEqual(met.annotations_source, {"test": ["manual", "model"]})

        # test the different modes to set sources
        met3 = MeMoMetabolite()
        # test names with lists
        met3.set_names(["Peter", "Hans"], source=["model", "manual"])
        self.assertEqual(met3.names_source, ["manual", "model"])
        # test annotations with dictionaries
        met3.set_annotations(met.annotations, met.annotations_source)
        self.assertEqual(met3.annotations_source, {"test": ["manual", "model"]})

        # test adding metabolites and ignoring duplicated metabolites
        met3.add_names(["Peter"], source=["duplicate"])
        self.assertEqual(met3.names_source, ["manual", "model"])
        met3.add_names(["Dieter"], source="no_dup")
        self.assertEqual(met3.names_source, ["no_dup", "manual", "model"])
        met3.add_annotations({"test": ["test"]}, source="duplicated")
        self.assertEqual(met3.annotations_source, {"test": ["manual", "model"]})
        met3.add_annotations({"test2": ["test"]}, source="no_dup")
        self.assertEqual(
            met3.annotations_source, {"test": ["manual", "model"], "test2": ["no_dup"]}
        )
        met3.add_inchi_string(
            "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1",
            source="pubchem",
        )
        self.assertEqual(met3._inchi_source, "pubchem")
        met3.add_inchi_string("test", source="duplicated")
        self.assertEqual(met3._inchi_source, "pubchem")

        # test some fail cases
        met = MeMoMetabolite()
        with self.assertRaises(ValueError):
            met.set_names(["Peter"], source=["model", "manual"])
        with self.assertRaises(ValueError):
            met.set_annotations({"A": ["a"]}, source={"B": ["model"]})
        with self.assertRaises(ValueError):
            met.set_annotations({"A": ["a"]}, source={"A": ["model", "test"]})


if __name__ == "__main__":
    unittest.main()
