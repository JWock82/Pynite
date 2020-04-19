from numpy import where, append, array, empty, float, delete, squeeze, size


class Debug():
    def debug(self):
        self.appNodes = array([['name', '1', '1', '1'],
                               ['name2', '2', '1', '1'],
                               ['name3', '3', '1', '1']])
        self.appMembers = array([['a', '0', '1', '0'], ['ab', '0', '2', '0']])
        self.appMaterials = array([['mat1', '1', '1']])

        mylist = array(['0', '2', '0'])
        a = self.duplicate(self.appMembers, mylist)
        print(a)

    def duplicate(self, array, find):
        r = False
        a = size(array, 1)
        i = 0
        while i <= len(array) - 1:
            try:
                j = 1
                while j <= a:
                    if (array[i, j] == find[j - 1]):
                        j += 1
                    else:
                        break
            except IndexError:
                if j == a:
                    r = True
            i += 1
            return r


ui = Debug()
ui.debug()
