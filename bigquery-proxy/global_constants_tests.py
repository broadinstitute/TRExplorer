import unittest

from global_constants import MAIN_BIGQUERY_TABLE_COLUMNS, find_problems_in_column_descriptions_shown_on_website


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


if __name__ == "__main__":
    unittest.main()
