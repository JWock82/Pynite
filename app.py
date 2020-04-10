# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 17:48:28 2020

@author: Tarang J. AUS
"""

from PyNite import FEModel3D
import tkinter as tk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import numpy as np
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow
import sys


def window():
    app = QApplication(sys.argv)
    win = QMainWindow()
    win.setGeometry(200, 200, 800, 400)
    win.setWindowTitle("PyNite Structural Analysis")

    win.show()
    sys.exit(app.exec_())


window()


class Main:
    def __init__(self, root, *args, **kwargs):

    def frameNode(self):
        self.nodeFrame = tk.Frame(master=root)
        self.nodeFrame.grid(row=0, column=9, columnspan=9)


class MainApplication:
    def __init__(self, root, truss, *args, **kwargs):
        self.root = root
        self.truss = truss
        self.main = Main(self.root)
        self.main.frameNode()
        self.nodes()
        self.members()

        # Application Storage
        self.appNodes = np.array(
            [['0987654321qwertyuiopasdfghjklzxcvbnm', 'x', 'y', 'z']],
            dtype="S")  # A list of the structure's nodes
        self.appauxNodes = []  # A list of the structure's auxiliary nodes
        self.appMembers = []  # A list of the structure's members
        self.appPlates = []  # A list of the structure's plates
        self.app__D = {
        }  # A dictionary of the structure's nodal displacements by load combination
        self.appLoadCombos = {
        }  # A dictionary of the structure's load combinations

    def nodeSubmit(self, name, x, y, z, oldname=None):
        name = str(name.strip())
        x = eval(x)
        y = eval(y)
        z = eval(z)
        if (oldname == None):
            oldname = name
        else:
            oldname.strip()

        self.nodeExists(oldname)

        if (self.exists == True):
            print("found")
            self.truss.RemoveNode(oldname)
            self.truss.AddNode(name, x, y, z)
            self.appNodes[self.check, 0] = name
            self.appNodes[self.check, 1] = str(x)
            self.appNodes[self.check, 2] = str(y)
            self.appNodes[self.check, 3] = str(z)
        if (self.exists == False):
            print("not found")
            self.truss.AddNode(name, x, y, z)
            self.appNodes = np.vstack(
                (self.appNodes, np.array([[name, str(x),
                                           str(y),
                                           str(z)]])))

        print(np.array([self.appNodes]))
        self.nodeDisplay()

    def nodeExists(self, name):
        exists = False
        check = np.squeeze(np.where(self.appNodes == name)[0])
        try:
            int(check)
        except TypeError:
            exists = False
        else:
            exists = True

        self.exists = exists
        self.check = check

    def nodeDelete(self, name):
        print(name)
        self.nodeExists(name)
        print(self.check)
        self.truss.RemoveNode(name)
        self.appNodes = np.delete(self.appNodes, self.check, 0)
        print(f"deleted {name}")
        self.nodeDisplay()

    def nodeDisplay(self):
        nodeTitle = tk.Label(self.main.nodeFrame, text="Nodes")
        nodeTitle.grid(row=0, column=9, columnspan=9)

        for i in range(len(self.appNodes) - 1):
            print(i)
            self.listnodeName = tk.Entry(self.main.nodeFrame)
            self.listnodeName.insert("end", self.appNodes[i + 1, 0])
            self.listnodeX = tk.Entry(self.main.nodeFrame)
            self.listnodeX.insert(0, self.appNodes[i + 1, 1])
            self.listnodeY = tk.Entry(self.main.nodeFrame)
            self.listnodeY.insert(0, self.appNodes[i + 1, 2])
            self.listnodeZ = tk.Entry(self.main.nodeFrame)
            self.listnodeZ.insert(0, self.appNodes[i + 1, 3])
            listnodeupdate = tk.Button(
                self.main.nodeFrame,
                text='Update',
                command=lambda: self.nodeSubmit(
                    self.listnodeName.get(),  # new name
                    self.listnodeX.get(),  # x
                    self.listnodeY.get(),  # y
                    self.listnodeZ.get(),  # z
                    self.appNodes[i + 1, 0]  # old name
                ))
            listnodedelete = tk.Button(
                self.main.nodeFrame,
                text='Delete',
                command=lambda: self.nodeDelete(self.appNodes[i + 1, 0]
                                                ))  # send name

            self.listnodeName.grid(row=i + 1, column=9)
            self.listnodeX.grid(row=i + 1, column=10)
            self.listnodeY.grid(row=i + 1, column=11)
            self.listnodeZ.grid(row=i + 1, column=12)
            listnodeupdate.grid(row=i + 1, column=13)
            listnodedelete.grid(row=i + 1, column=14)
        i = 0

    def updateFrame(self):
        # Analyze the model
        self.truss.Analyze()

    def nodes(self):
        # Add Nodes
        nodeLabel = tk.Label(root, text="Add Node:")
        nodeLabel.grid(row=1, column=0)

        nodeName = tk.Entry(root)
        nodeName.insert(0, "Node name")
        nodeX = tk.Entry(root)
        nodeX.insert(0, "x (m)")
        nodeY = tk.Entry(root)
        nodeY.insert(0, "y (m)")
        nodeZ = tk.Entry(root)
        nodeZ.insert(0, "z (m)")

        nodeButton = tk.Button(
            root,
            text="Add Node",
            command=lambda: self.nodeSubmit(nodeName.get(), nodeX.get(),
                                            nodeY.get(), nodeZ.get()))

        nodeName.grid(row=1, column=1)
        nodeX.grid(row=1, column=2)
        nodeY.grid(row=1, column=3)
        nodeZ.grid(row=1, column=4)
        nodeButton.grid(row=1, column=5)

    # Add Supports

    def members(self):
        # Add Members

        memberLabel = tk.Label(self.root, text="Create Member:")
        memberLabel.grid(row=2, column=0)

        memberName = tk.Entry(self.root)
        memberName.insert(0, "Member name")

        self.memberIVar = tk.StringVar(root)
        self.memberIVar.set("Select I node")
        self.memberJVar = tk.StringVar(root)
        self.memberJVar.set("Select J node")

        nodelist = [""]
        self.memberI = tk.OptionMenu(self.root, self.memberIVar, *nodelist)
        self.memberJ = tk.OptionMenu(self.root, self.memberJVar, *nodelist)

        memberButton = tk.Button(root, text="Add Member")

        memberName.grid(row=2, column=1)
        self.memberI.grid(row=2, column=2)
        self.memberJ.grid(row=2, column=3)
        memberButton.grid(row=2, column=4)

    def memberListUpdate(self):
        nodelist = []

        for i in range(len(self.appNodes) - 1):
            nodelist.append(self.appNodes[i + 1, 0])
        print(nodelist)
        self.memberI = tk.OptionMenu(self.root, self.memberIVar, *nodelist)
        self.memberJ = tk.OptionMenu(self.root, self.memberJVar, *nodelist)


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("1280x720")
    truss = FEModel3D()
    MainApplication(root, truss)
    root.mainloop()

# Add Forces

# Create members (all members will have the same properties in this example)
J = 250
Iy = 250
Iz = 200
E = 30000
G = 250
A = 12

# frame.AddMember('M1', 'N1', 'N2', E, G, Iy, Iz, J, A)
# frame.AddMember('M2', 'N2', 'N3', E, G, Iy, Iz, J, A)
# frame.AddMember('M3', 'N3', 'N4', E, G, Iy, Iz, J, A)
# frame.AddMember('M4', 'N4', 'N5', E, G, Iy, Iz, J, A)
# frame.AddMember('M5', 'N5', 'N6', E, G, Iy, Iz, J, A)

# Add nodal loads
# frame.AddNodeLoad('N3', 'FY', -30)
# frame.AddNodeLoad('N4', 'FY', -30)

# node1 = frame.GetNode('N1')
# node6 = frame.GetNode('N6')

# print('Calculated reactions: ', node1.RxnFX, node1.RxnFY, node1.RxnMZ,
#       node6.RxnFX, node6.RxnFY, node6.RxnMZ)
# print('Expected reactions: ', 11.69, 30, -1810, -11.69, 30, 1810)
# print('Calculated displacements: ',
#       frame.GetNode('N3').DY,
#       frame.GetNode('N4').DY,
#       frame.GetNode('N3').RZ,
#       frame.GetNode('N4').RZ)
