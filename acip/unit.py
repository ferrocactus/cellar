class Unit:
    def __init__(self, verbose=False, **args):
        self._verbose = verbose
        self.args = args

    def get(self):
        pass

    def vprint(self, *args):
        if self._verbose == True:
            print(*args)