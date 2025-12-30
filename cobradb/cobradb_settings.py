# -*- coding: utf-8 -*-

"""Retrieve local user settings"""

from configparser import ConfigParser, NoOptionError
import os
from os.path import normpath, isfile, expanduser

config = ConfigParser()


class CobraDBSettings:
    postgres_host: str = "localhost"
    postgres_port: str = "5432"
    postgres_user: str = "cobradb"
    postgres_password: str = ""
    postgres_database: str = "cobradb"
    postgres_test_database: str = "cobradb_test"

    compartment_names: str = ""
    data_source_preferences: str = ""
    rhea: str = ""
    curated_metabolites: str = ""
    curated_reactions: str = ""
    taxonomy_nodes: str = ""
    taxonomy_names: str = ""

    n_processes: int = 1

    def __init__(self):
        # overwrite defaults settings with settings from the file
        # filepath = abspath(join(dirname(__file__), '..', 'settings.ini'))
        filepath = normpath("settings.ini")
        if isfile(filepath):
            config.read(filepath)
        else:
            raise Exception("No settings files at path: %s" % filepath)

        # prefer environment variables for database settings
        env_names = {
            "postgres_host": "COBRADB_POSTGRES_HOST",
            "postgres_port": "COBRADB_POSTGRES_PORT",
            "postgres_user": "COBRADB_POSTGRES_USER",
            "postgres_password": "COBRADB_POSTGRES_PASSWORD",
            "postgres_database": "COBRADB_POSTGRES_DATABASE",
            "postgres_test_database": "COBRADB_POSTGRES_TEST_DATABASE",
        }
        for setting_name, env_name in env_names.items():
            if env_name in os.environ:
                print(
                    "Setting %s with environment variable %s" % (setting_name, env_name)
                )
                setattr(self, setting_name, os.environ[env_name])
            else:
                setattr(self, setting_name, config.get("DATABASE", setting_name))

        # set up the database connection string
        self.db_connection_string = "postgresql://%s:%s@%s:%s/%s" % (
            self.postgres_user,
            self.postgres_password,
            self.postgres_host,
            self.postgres_port,
            self.postgres_database,
        )

        if config.has_option("GENERAL", "n_processes"):
            self.n_processes = int(config.get("GENERAL", "n_processes"))

        if not config.has_section("DATA"):
            raise Exception("DATA section was not found in settings.ini")

        # these are required
        try:
            self.model_directory = expanduser(config.get("DATA", "model_directory"))
        except NoOptionError:
            raise Exception("model_directory was not supplied in settings.ini")

        try:
            self.refseq_directory = expanduser(config.get("DATA", "refseq_directory"))
        except NoOptionError:
            raise Exception("refseq_directory was not supplied in settings.ini")
        try:
            self.model_genome = expanduser(config.get("DATA", "model_genome"))
        except NoOptionError:
            raise Exception("model_genome path was not supplied in settings.ini")

        # these are optional
        for data_pref in [
            "compartment_names",
            "data_source_preferences",
            "rhea",
            "curated_metabolites",
            "curated_reactions",
            "taxonomy_nodes",
            "taxonomy_names",
            "chebi_cache_dir",
            "default_chebis",
        ]:
            try:
                setattr(self, data_pref, expanduser(config.get("DATA", data_pref)))
            except NoOptionError:
                setattr(self, data_pref, None)


settings = CobraDBSettings()
