import unittest

from global_constants import MAIN_BIGQUERY_TABLE_COLUMNS, GROUP_ROWS_IN_COLUMN_DIALOGS, PLOT_AXIS_TYPE_CATEGORICAL, \
    PLOT_AXIS_TYPE_CATEGORICAL_OR_NUMERICAL, PLOT_AXIS_TYPE_NUMERICAL, find_problems_in_column_descriptions_shown_on_website, \
    get_plot_columns


def column(description):
    return {"name": "TestColumn", "description": description}


class FindProblemsInColumnDescriptionsShownOnWebsiteTests(unittest.TestCase):

    def test_allowed_tags_and_a_less_than_sign_followed_by_a_space_pass(self):
        description = ("Values:<ul><li><b>“significant”</b>: FDR < 0.1</li></ul><br /><i>see</i> "
                       "<a href='https://example.org'>source</a>, <code>1:100-200</code>")
        self.assertEqual(find_problems_in_column_descriptions_shown_on_website([column(description)]), [])

    def test_literal_angle_bracket_text_is_reported_as_an_unknown_tag(self):
        problems = find_problems_in_column_descriptions_shown_on_website(
            [column("Inner span of the variation cluster (<VC:...>), or empty for an isolated TR (<TR:...>)")])
        self.assertEqual(len(problems), 2)
        self.assertIn("<tr:...>", problems[0])
        self.assertIn("<vc:...>", problems[1])

    def test_escaped_angle_brackets_pass(self):
        self.assertEqual(find_problems_in_column_descriptions_shown_on_website(
            [column("Inner span of the variation cluster (&lt;VC:...&gt;)")]), [])

    def test_straight_double_quote_is_reported(self):
        problems = find_problems_in_column_descriptions_shown_on_website([column('The "All" view')])
        self.assertEqual(problems, ["TestColumn: contains a straight double quote; use “ and ” instead"])

    def test_column_without_a_description_passes(self):
        self.assertEqual(find_problems_in_column_descriptions_shown_on_website([{"name": "TestColumn"}]), [])

    def test_main_table_column_descriptions_have_no_problems(self):
        self.assertEqual(find_problems_in_column_descriptions_shown_on_website(MAIN_BIGQUERY_TABLE_COLUMNS), [])


class GroupRowsInColumnDialogsTests(unittest.TestCase):

    def test_every_column_group_appears_in_exactly_one_dialog_row(self):
        """A group missing from the rows would leave its columns out of the export and custom filter dialogs."""
        groups_in_rows = [group for row in GROUP_ROWS_IN_COLUMN_DIALOGS for group in row]
        self.assertEqual(len(groups_in_rows), len(set(groups_in_rows)))
        self.assertEqual({c["group"] for c in MAIN_BIGQUERY_TABLE_COLUMNS if c.get("group")}, set(groups_in_rows))


class GetPlotColumnsTests(unittest.TestCase):

    def test_axis_type_defaults_to_the_one_for_the_column_type(self):
        columns = [
            {"name": "S", "type": "STRING", "allowCustomFilter": True},
            {"name": "I", "type": "INTEGER", "allowCustomFilter": True},
            {"name": "F", "type": "FLOAT", "allowCustomFilter": True},
        ]
        self.assertEqual([c["plotAxisType"] for c in get_plot_columns(columns)],
                         [PLOT_AXIS_TYPE_CATEGORICAL, PLOT_AXIS_TYPE_CATEGORICAL_OR_NUMERICAL, PLOT_AXIS_TYPE_NUMERICAL])

    def test_plot_axis_type_key_overrides_the_default_and_none_leaves_the_column_out(self):
        columns = [
            {"name": "Start", "type": "INTEGER", "allowCustomFilter": True, "plotAxisType": PLOT_AXIS_TYPE_NUMERICAL},
            {"name": "LocusId", "type": "STRING", "allowCustomFilter": True, "plotAxisType": None},
        ]
        self.assertEqual(get_plot_columns(columns), [
            {"name": "Start", "displayName": "Start", "group": "", "type": "INTEGER",
             "plotAxisType": PLOT_AXIS_TYPE_NUMERICAL},
        ])

    def test_columns_without_custom_filters_are_left_out(self):
        self.assertEqual(get_plot_columns([{"name": "id", "type": "INTEGER"}]), [])

    def test_main_table_identifiers_are_left_out_and_motif_size_can_be_either(self):
        axis_types = {c["name"]: c["plotAxisType"] for c in get_plot_columns()}
        self.assertNotIn("LocusId", axis_types)
        self.assertNotIn("ReferenceRepeatSequence", axis_types)
        self.assertEqual(axis_types["MotifSize"], PLOT_AXIS_TYPE_CATEGORICAL_OR_NUMERICAL)
        self.assertEqual(axis_types["GencodeGeneRegion"], PLOT_AXIS_TYPE_CATEGORICAL)
        self.assertEqual(axis_types["HPRC256_Stdev"], PLOT_AXIS_TYPE_NUMERICAL)


if __name__ == "__main__":
    unittest.main()
