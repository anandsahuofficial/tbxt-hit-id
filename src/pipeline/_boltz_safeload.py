"""
Boltz launcher that monkey-patches torch.load to weights_only=False before
boltz tries to load its checkpoint. Required for torch 2.6+ where the default
flipped to weights_only=True and breaks Boltz's omegaconf-pickled checkpoints.
"""
import sys

import torch
_orig_load = torch.load

def _patched_load(*args, **kwargs):
    kwargs["weights_only"] = False
    return _orig_load(*args, **kwargs)

torch.load = _patched_load

try:
    import torch.serialization as _ts
    from omegaconf.dictconfig import DictConfig
    from omegaconf.listconfig import ListConfig
    from omegaconf.base import ContainerMetadata, Metadata
    _ts.add_safe_globals([DictConfig, ListConfig, ContainerMetadata, Metadata])
except Exception:
    pass

from boltz.main import cli

if __name__ == "__main__":
    sys.argv[0] = "boltz"
    cli()
