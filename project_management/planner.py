from project_management.configuration import Configuration

class Planner():
    
    @staticmethod
    def run(args: list[str], config: Configuration) -> None:
        
        # Commencing pipeline
        if args.fresh:
            # TODO clear manifests and data directories
            pass
        
        # TODO figure out where to begin pipeline and start there
            
    #==========================================================================
    # Private Methods Below
    #==========================================================================
        
    @staticmethod
    def _