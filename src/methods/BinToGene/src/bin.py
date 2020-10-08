class Bin:
    """
    A single bin. Contains sequence name, start, end, index.
    """
    def __init__(self, seqname, start, end, index):
        """
        Parameters:
        ___________

        seqname: sequence name, i.e., chromosome

        start: start location

        end: end location

        index: index in the original cell-by-bin matrix
        """
        self.seqname = seqname
        if start >= end:
            raise ValueError(f"Invalid bin interval ({start}, {end})")
        self.start = start
        self.end = end
        self.index = index

    def intersects(self, start, end):
        """
        Check whether bin intersects the given (start, end) range.
        """
        return not ((end < self.start) or (self.end < start))

    def precedes(self, start, end):
        """
        Check whether bin precedes the given (start, end) range.
        """
        return self.end < start

    def suceeds(self, start, end):
        """
        Check whether bin suceeds the given (start, end) range.
        """
        return end < self.start

