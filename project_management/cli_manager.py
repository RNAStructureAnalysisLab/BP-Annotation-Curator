# AUTHOR: Kristopher Church

import argparse

'''
PURPOSE: Parse command line arguments for the pipeline
USE: CLIManager.run() returns parsed arguments
'''
class CLIManager:
    DEFAULT_VERSION = "3.95"
    
    @staticmethod
    def run(args: list[str]=None) -> argparse.Namespace:
        parser = CLIManager._build_parser()
        parsed = parser.parse_args(args)
        
        return parsed
        
    #==========================================================================
    # HELPER METHODS BELOW
    #==========================================================================
    
    @staticmethod
    def _build_parser() -> argparse.ArgumentParser:
        parser = argparse.ArgumentParser(
            description="RNA Base Pairing Annotation Pipeline",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog=
                '''
                Examples:
                    python main.py                          Run pipeline (default using version 3.95)
                    python main.py -v 3.95                  Run pipeline for version 3.95
                    python main.py --fresh          Start fresh, ignore existing data
                    python main.py --dry-run        Show what would run without executing
                '''
        )
            
        parser.add_argument(
            "-v", "--version",
            type=str,
            default="3.95",
            help=f"Data version to download/process (default: {CLIManager.DEFAULT_VERSION})"
        )
        
        parser.add_argument(
            "--fresh",
            action="store_true",
            default=False,
            help="Start pipeline from beginning, ignoring existing data"
        )
        
        parser.add_argument(
            "--dry-run",
            action="store_true",
            default=False,
            help="Show what steps would run without executing them"
        )
        
        return parser