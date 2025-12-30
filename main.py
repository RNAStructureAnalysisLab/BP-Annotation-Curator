# AUTHOR: Kristopher Church

import os
import sys

from project_management.configuration import Configuration
from project_management.cli_manager import CLIManager
from project_management.planner import Planner
from project_management.manifest_manager import ManifestManager
from data_management.raw_data_loading.cluster_downloader import ClusterDownloader
from data_management.raw_data_loading.pdb_downloader import PdbDownloader

def main() -> None:
    parsed_args = CLIManager.run(sys.argv)
    config = Configuration(parsed_args.version)
    Planner.run(parsed_args, config)

if __name__ == "__main__":
    main()
    