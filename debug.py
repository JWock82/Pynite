import numpy as np
import tkinter as tk


class Debug():
    def __init__(self):
        self.truss = np.array([['dont show this row', 'a', 'a', 'a']],
                              dtype="S")
        self.truss = np.vstack(
            (self.truss, np.array([['Node name', '1', '2', '3']])))
        self.truss = np.vstack(
            (self.truss, np.array([['Node name 1', '1', '2', '3']])))

    def debug(self):
        print(self.truss)
        # oldname = "Node name"
        # newname = "Node name"
        # x = str(41)
        # y = str(43)
        # z = str(4)
        # exists = False
        # check = np.squeeze(np.where(self.truss == oldname)[0])
        # try:
        #     int(check)
        # except TypeError:
        #     exists = False
        # else:
        #     exists = True

        # if (exists == True):
        #     print(check)
        #     print("is true")
        #     self.truss[check, 0] = newname
        #     self.truss[check, 1] = x
        #     self.truss[check, 2] = y
        #     self.truss[check, 3] = z

        # else:
        #     print("not found")
        #     self.truss = np.vstack(
        #         (self.truss, np.array([[newname,
        #                                 str(x),
        #                                 str(y),
        #                                 str(z)]])))
        self.truss = np.delete(self.truss, 1, 0)
        print(self.truss)


start = Debug()
start.debug()
