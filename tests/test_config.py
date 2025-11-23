from pathlib import Path

import pandas as pd
import pytest

import staging_visualization as sv


class TestDeepMergeDict:
    def test_shallow_merge(self):
        base = {"a": 1, "b": 2}
        override = {"b": 3, "c": 4}

        result = sv.deep_merge_dict(base, override)

        assert result == {"a": 1, "b": 3, "c": 4}

    def test_nested_merge(self):
        base = {"outer": {"inner": 1, "keep": 2}}
        override = {"outer": {"inner": 99}}

        result = sv.deep_merge_dict(base, override)

        assert result["outer"]["inner"] == 99
        assert result["outer"]["keep"] == 2

    def test_base_not_mutated(self):
        base = {"a": {"b": 1}}
        override = {"a": {"b": 2}}

        sv.deep_merge_dict(base, override)

        assert base["a"]["b"] == 1


class TestLoadConfig:
    def test_missing_config_raises_systemexit(self):
        with pytest.raises(SystemExit):
            sv.load_config(Path("/nonexistent/config.yaml"))

    def test_valid_config_returns_expected_keys(self, tmp_path):
        config_file = tmp_path / "config.yaml"
        config_file.write_text(
            """
file_settings:
  default_data_path: data/test.tsv
  default_output_dir: output/
  figure_dpi: 150
""",
            encoding="utf-8",
        )

        config = sv.load_config(config_file)

        assert "file_settings" in config
        assert "column_mapping" in config
        assert "visualization" in config

    def test_invalid_yaml_exits(self, tmp_path, capsys):
        config_file = tmp_path / "config.yaml"
        config_file.write_text(":\n- bad", encoding="utf-8")

        with pytest.raises(SystemExit):
            sv.load_config(config_file)

        captured = capsys.readouterr()
        assert "Failed to parse configuration file" in captured.err

    def test_non_mapping_config_exits(self, tmp_path, capsys):
        config_file = tmp_path / "config.yaml"
        config_file.write_text("- just\n- a\n- list\n", encoding="utf-8")

        with pytest.raises(SystemExit):
            sv.load_config(config_file)

        assert "must contain a mapping" in capsys.readouterr().err

    @pytest.mark.parametrize(
        "key,content,expected_message",
        [
            ("file_settings", "file_settings: invalid", "file_settings must be a mapping"),
            ("column_mapping", "column_mapping: invalid", "column_mapping must be a mapping"),
            ("required_columns", "required_columns: invalid", "required_columns must be a list"),
            ("stage_ordering", "stage_ordering: invalid", "stage_ordering must be a mapping"),
            ("visualization", "visualization: invalid", "visualization must be a mapping"),
        ],
    )
    def test_invalid_sections_exit(self, tmp_path, capsys, key, content, expected_message):
        base_config = """
file_settings:
  default_data_path: data/test.tsv
  default_output_dir: output/
  figure_dpi: 150
"""
        config_file = tmp_path / f"{key}.yaml"
        config_file.write_text(base_config + content + "\n", encoding="utf-8")

        with pytest.raises(SystemExit):
            sv.load_config(config_file)

        assert expected_message in capsys.readouterr().err


class TestResolveAndLoadCohort:
    def test_resolve_path_relative_and_absolute(self, tmp_path):
        absolute = sv.resolve_path(tmp_path, relative_to_base=False)
        assert absolute == tmp_path

        relative = sv.resolve_path("relative/path.txt", relative_to_base=True)
        assert relative.is_absolute()

    def test_load_cohort_missing_file_exits(self, tmp_path, capsys):
        missing = tmp_path / "missing.tsv"

        with pytest.raises(SystemExit):
            sv.load_cohort(missing)

        assert "Could not locate cohort file" in capsys.readouterr().err

    def test_load_cohort_parser_error_exits(self, tmp_path, capsys, monkeypatch):
        bad_file = tmp_path / "bad.tsv"
        bad_file.write_text("a\tb\n1\t2\t3\n", encoding="utf-8")

        def fake_read_csv(*args, **kwargs):
            raise pd.errors.ParserError("boom")

        monkeypatch.setattr(pd, "read_csv", fake_read_csv)

        with pytest.raises(SystemExit):
            sv.load_cohort(bad_file)

        assert "Could not parse cohort file" in capsys.readouterr().err
