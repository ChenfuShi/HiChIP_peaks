from context import hichip_tool
#the error above is fine when i run it from vscode and command line
#probably can't realize the path modification

from hichip_tool import interaction_to_sparse

import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s - %(message)s",
    handlers=[
    logging.StreamHandler()
]
)
logging.error("ciao")