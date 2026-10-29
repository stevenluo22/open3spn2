class DNAForce(object):
    """ Wrapper for the openMM force. """
    _table = None

    def __init__(self, dna):
        self.periodic = dna.periodic
        self.force = None
        self.dna = dna

        # Define the dna force
        self.reset()

        # Define the interaction pairs
        self.defineInteraction()

    def __getattr__(self, attr):
        if attr in self.__dict__:
            return getattr(self, attr)
        elif 'force' in self.__dict__:
            return getattr(self.force, attr)
        else:
            if '__repr__' in self.__dict__:
                raise AttributeError(f"type object {str(self)} has no attribute {str(attr)}")
            else:
                raise AttributeError()

    def addForce(self, system):
        """Add this force's OpenMM force to `system`. Forces that manage several OpenMM forces
        (e.g. BasePair, CrossStacking) override this to add all of them.

        The system takes ownership of the OpenMM force, so the interaction table is read before it is added."""
        if hasattr(type(self), 'read_table'):
            self.table
        system.addForce(self.force)

    @property
    def table(self):
        """DataFrame of the interactions of this force (particle indices and parameters), read once from the OpenMM
        force(s) by `read_table` and cached."""
        if self._table is None:
            self._table = self.read_table()
        return self._table

    def _checked_box(self, box):
        if self.periodic and box is None:
            raise ValueError(f'{type(self).__name__} uses periodic boundary conditions: a box is required')
        if not self.periodic and box is not None:
            raise ValueError(f'{type(self).__name__} does not use periodic boundary conditions: no box expected')
        return box


class ProteinDNAForce(DNAForce):
    def __init__(self, dna, protein):
        self.protein = protein
        super().__init__(dna)

