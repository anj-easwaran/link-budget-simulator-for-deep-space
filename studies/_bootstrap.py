import os
import sys

# full path to the model folder, which sits one level up from this studies folder
_MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          os.pardir, "model")

# adds model/ to the front of the search path, but only if it isn't there already
if _MODEL_DIR not in sys.path:
    sys.path.insert(0, _MODEL_DIR)
