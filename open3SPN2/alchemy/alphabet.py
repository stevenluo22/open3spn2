"""The names a base-pair position of a duplex can take.

An identity is the pair of bases (strand A, strand B) of a position. A mix is a name that is not an identity: it
has the parameters of identity a + weight (b - a), interpolated jointly at every position of an interaction where the
mix appears. Names are ordered identities first, then mixes; a state is an array with the index of the name of each
position.
"""
WATSON_CRICK = {'GC': ('G', 'C'), 'AT': ('A', 'T'), 'TA': ('T', 'A'), 'CG': ('C', 'G')}


class Alphabet:
    """Identities {name: (base on strand A, base on strand B)} and mixes {name: (identity a, identity b, weight)}."""

    def __init__(self, identities, mixes=None):
        self.identities = {name: tuple(pair) for name, pair in identities.items()}
        self.mixes = {name: (a, b, float(weight)) for name, (a, b, weight) in (mixes or {}).items()}
        for name, (a, b, weight) in self.mixes.items():
            if name in self.identities or a not in self.identities or b not in self.identities or a == b:
                raise ValueError(f'mix {name} must interpolate between two identities and not reuse an identity name')
            if not 0 < weight < 1:
                raise ValueError(f'mix {name} has weight {weight}, outside (0, 1)')
        self.names = list(self.identities) + list(self.mixes)

    def __len__(self):
        return len(self.names)

    def index(self, name):
        return self.names.index(name)

    def __repr__(self):
        return f'Alphabet({self.names})'


def alphabet(preset='wc'):
    """The alphabet of a preset: 'wc' (GC, AT, TA, CG and W halfway between AT and TA) or 'wc4' (the four
    Watson-Crick identities)."""
    if preset == 'wc':
        return Alphabet(WATSON_CRICK, {'W': ('AT', 'TA', 0.5)})
    if preset == 'wc4':
        return Alphabet(WATSON_CRICK)
    raise ValueError(f'unknown alphabet preset {preset!r}')
