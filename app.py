# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 17:48:28 2020

@author: Tarang J. AUS =)
"""

# Required PyNite Module
from PyNite import FEModel3D
# PyQt5 used for GUI interface - built in QT Designer Software
from PyQt5 import QtCore, QtGui, QtWidgets
# Get .ui file (converted to .py) for GUI
from design import Ui_MainWindow
# One use-case of 'lambda' does not work, so using partial
from functools import partial
# For general math functions
from numpy import (empty, vstack, array, delete, squeeze, where, float)
# For application
from sys import argv
# For plotting
import matplotlib
matplotlib.use('qt5agg')
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg,
                                                NavigationToolbar2QT as
                                                NavigationToolbar)
from matplotlib.figure import Figure
# For special functions as inputs
func = '''(Abs, acos, acosh, acot, acoth, acsc, acsch, adjoint, airyai,
            airyaiprime, airybi, airybiprime, appellf1, arg, asec, asech, asin,
            asinh, assoc_laguerre, assoc_legendre, atan, atan2, atanh, bell,
            bernoulli, besseli, besselj, besselk, bessely, beta, binomial,
            bspline_basis, bspline_basis_set, carmichael, catalan, cbrt,
            ceiling, chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
            Chi, Ci, combinatorial, conjugate, cos, cosh, cot, coth, csc, csch,
            digamma, DiracDelta, dirichlet_eta, E1, Ei, Eijk, elementary,
            elliptic_e, elliptic_f, elliptic_k, elliptic_pi, erf, erf2,
            erf2inv, erfc, erfcinv, erfi, erfinv, euler, exp, exp_polar,
            expint, factorial, factorial2, FallingFactorial, ff, fibonacci,
            floor, frac, fresnelc, fresnels, gamma, gegenbauer, genocchi,
            hankel1, hankel2, harmonic, Heaviside, hermite, hn1, hn2, hyper,
            Id, im, interpolating_spline, jacobi, jacobi_normalized, jn,
            jn_zeros, KroneckerDelta, laguerre, LambertW, legendre, lerchphi,
            LeviCivita, Li, li, ln, log, loggamma, lowergamma, lucas, marcumq,
            mathieuc, mathieucprime, mathieus, mathieusprime, Max, meijerg,
            Min, multigamma, partition, periodic_argument, Piecewise,
            piecewise_fold, polar_lift, polarify, polygamma, polylog,
            principal_branch, re, real_root, rf, RisingFactorial, root, sec,
            sech, Shi, Si, sign, sin, sinc, SingularityFunction, sinh, special,
            sqrt, stieltjes, subfactorial, tan, tanh, transpose, tribonacci,
            trigamma, unbranched_argument, unpolarify, uppergamma, yn, Ynm,
            Ynm_c, zeta, Znm)'''
from sympy import E, pi, I
exec("from sympy.functions import " + func)


class ApplicationWindow(Ui_MainWindow):
    def __init__(self, MainWindow, truss, *args, **kwargs):

        # Variable used to communicate with FEModel3D PyNite Program.
        self.truss = truss

        #---------------------------------------------------------------------#
        #                        Initiate GUI Interface                       #
        #---------------------------------------------------------------------#

        self.MainWindow = MainWindow
        self.design = Ui_MainWindow(
        )  # use generated design file from PyQt Designer
        self.design.setupUi(self.MainWindow)

        # Member Plot
        # memberFig = Figure(dpi=dpi)
        # self.member_ax = memberFig.add_subplot(projection='3d')

        #---------------------------------------------------------------------#
        #                        Application Storage                          #
        #---------------------------------------------------------------------#

        self.appNodes = empty((1, 4))  # A list of the structure's nodes
        self.appMembers = empty((1, 3))  # A list of the structure's members
        self.appPlates = []  # A list of the structure's plates
        self.app__D = {
        }  # A dictionary of the structure's nodal displacements by load combination
        self.appLoadCombos = {
        }  # A dictionary of the structure's load combinations

        # Variable to check if deleted empty row from 'self.appNodes'
        self.alreadyDeletedn = False
        self.alreadyDeletedm = False

        #---------------------------------------------------------------------#
        #                     Initiate Connect Functions                      #
        #---------------------------------------------------------------------#

        self.nodeStart()
        self.memberStart()

    def updateFEModel3D(self):
        # Analyze the model
        self.truss.Analyze()

    #-------------------------------------------------------------------------#
    #                                                                         #
    #                        Button Connect Functions                         #
    #                                                                         #
    #-------------------------------------------------------------------------#

    def nodeStart(self):
        # Link node button to 'self.nodeSubmit'
        # Inputs are taken after the button is pressed so that the latest values are used.
        self.design.addNodeButton.clicked.connect(
            lambda: self.nodeSubmit(self.design.defineNodeName.text().strip(
            ), self.design.defineXnode.text(), self.design.defineYnode.text(),
                                    self.design.defineZnode.text()))

        # Initiate plot canvas to allow updating values
        # self.nodePlot()

    def memberStart(self):
        # Link member button to 'self.memberSubmit'
        # Inputs are taken after the button is pressed so that the latest values are used.
        self.design.addMemberButton.clicked.connect(lambda: self.memberSubmit(
            self.design.defineMemberName.text().strip(),
            self.design.defineNodeI.text(), self.design.defineYnode.text(),
            self.design.defineZnode.text()))

        # Initiate plot canvas to allow updating values
        # self.memberPlot()

    #-------------------------------------------------------------------------#
    #                                                                         #
    #                             Common Functions                            #
    #                                                                         #
    #-------------------------------------------------------------------------#

    def tableClear(self, table):
        #---------------------------------------------------------------------#
        #                        Delete table contents                        #
        #---------------------------------------------------------------------#

        # Delete all contents inside table
        for i in reversed(range(table.count())):
            layout.itemAt(i).widget().setParent(None)

    #-------------------------------------------------------------------------#
    #                                                                         #
    #                             Node Functions                              #
    #                                                                         #
    #-------------------------------------------------------------------------#

    def nodeSubmit(self, name, x, y, z, oldname=""):
        # Check if input is valid
        if name != "" and x != "" and y != "" and z != "":
            try:
                x = float(eval(x))
                y = float(eval(y))
                z = float(eval(z))
            except NameError:
                print("Coordinate input must contain valid real numbers.")
                return
            except TypeError:
                print("Coordinate input must contain valid real numbers.")
                return
            except SyntaxError:
                print("Syntax Error.")
                return
            except RecursionError:
                print("Cannot contain recursive input.")
                return
            except (OverflowError, FloatingPointError):
                print("Overflow Error")
                return
            except ZeroDivisionError:
                print("Cannot divide by zero.")
                return
            else:
                print("Node definition passed.")
        else:
            print("Please enter a value.")
            return

        # The following is for when user edits the displayed text (mistakes)
        if (oldname == ""):
            oldname = name
        else:
            oldname.strip()

        self.nodeExists(oldname)

        print("Node name:", name)
        print("x:", x)
        print("y:", y)
        print("z:", z)

        # If exists, delete and replace; else, create new.
        if (self.exists == True):
            print("Node already exists.")
            self.truss.RemoveNode(oldname)
            self.truss.AddNode(name, x, y, z)
            self.appNodes[self.check, 0] = name
            self.appNodes[self.check, 1] = str(x)
            self.appNodes[self.check, 2] = str(y)
            self.appNodes[self.check, 3] = str(z)
        if (self.exists == False):
            print("Node does not exist.")
            self.truss.AddNode(name, x, y, z)
            self.appNodes = vstack(
                (self.appNodes, array([name, str(x),
                                       str(y), str(z)])))

        print("-" * 80)

        # Remove empty row in 'self.appNodes' (empty array)
        if (self.alreadyDeletedn == False):
            self.appNodes = delete(self.appNodes, 0, 0)
            self.alreadyDeletedn = True

        print(self.appNodes)
        self.nodeDisplay()

    def nodeExists(self, name):
        #---------------------------------------------------------------------#
        #                         Check if node exists                        #
        #---------------------------------------------------------------------#

        exists = False
        check = squeeze(where(self.appNodes == name)[0])
        try:
            int(check)
        except TypeError:
            exists = False
        else:
            exists = True

        self.exists = exists
        self.check = check

    def nodeDisplay(self):
        #---------------------------------------------------------------------#
        #                            Display nodes                            #
        #---------------------------------------------------------------------#

        i = 0
        for i in range(len(self.appNodes)):
            # set variables
            self.scrollAreaWidgetContents = self.design.scrollAreaWidgetContents
            self.table2editNode = self.design.table2editNode

            # node name column
            exec("self.editNodeName_" + str(i) +
                 " = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)")
            exec("self.editNodeName_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(200, 0))")
            exec("self.editNodeName_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(400, 16777215))")
            exec("self.editNodeName_" + str(i) +
                 ".setObjectName('editNodeName_" + str(i) + "')")
            exec("self.editNodeName_" + str(i) + ".setText('" +
                 self.appNodes[i, 0] + "')")
            exec("self.table2editNode.addWidget(self.editNodeName_" + str(i) +
                 ", " + str(i) + ", 0, 1, 1)")

            # node x column
            exec("self.editXnode_" + str(i) +
                 " = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)")
            exec("self.editXnode_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.editXnode_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(100, 16777215))")
            exec("self.editXnode_" + str(i) + ".setObjectName('editXnode_" +
                 str(i) + "')")
            exec("self.editXnode_" + str(i) + ".setText('" +
                 str(self.appNodes[i, 1]) + "')")
            exec("self.table2editNode.addWidget(self.editXnode_" + str(i) +
                 ", " + str(i) + ", 1, 1, 1)")

            # node y column
            exec("self.editYnode_" + str(i) +
                 " = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)")
            exec("self.editYnode_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.editYnode_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(100, 16777215))")
            exec("self.editYnode_" + str(i) + ".setObjectName('editYnode_" +
                 str(i) + "')")
            exec("self.editYnode_" + str(i) + ".setText('" +
                 self.appNodes[i, 2] + "')")
            exec("self.table2editNode.addWidget(self.editYnode_" + str(i) +
                 ", " + str(i) + ", 2, 1, 1)")

            # node z column
            exec("self.editZnode_" + str(i) +
                 " = QtWidgets.QLineEdit(self.scrollAreaWidgetContents)")
            exec("self.editZnode_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.editZnode_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(100, 16777215))")
            exec("self.editZnode_" + str(i) + ".setObjectName('editZnode_" +
                 str(i) + "')")
            exec("self.editZnode_" + str(i) + ".setText('" +
                 self.appNodes[i, 3] + "')")
            exec("self.table2editNode.addWidget(self.editZnode_" + str(i) +
                 ", " + str(i) + ", 3, 1, 1)")

            # node update button
            exec("self.updateNodeButton_" + str(i) +
                 " = QtWidgets.QPushButton(self.scrollAreaWidgetContents)")
            exec("self.updateNodeButton_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(90, 16777215))")
            exec("self.updateNodeButton_" + str(i) +
                 ".setObjectName('updateNodeButton_" + str(i) + "')")
            exec("self.table2editNode.addWidget(self.updateNodeButton_" +
                 str(i) + ", " + str(i) + ", 4, 1, 1)")

            # node delete button
            exec("self.deleteNodeButton_" + str(i) +
                 " = QtWidgets.QPushButton(self.scrollAreaWidgetContents)")
            exec("self.deleteNodeButton_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(90, 16777215))")
            exec("self.deleteNodeButton_" + str(i) +
                 ".setObjectName('deleteNodeButton_" + str(i) + "')")
            exec("self.table2editNode.addWidget(self.deleteNodeButton_" +
                 str(i) + ", " + str(i) + ", 5, 1, 1)")

            exec(
                "self.updateNodeButton_" + str(i) +
                ".setText(QtCore.QCoreApplication.translate('MainWindow','Update'))"
            )
            exec(
                "self.deleteNodeButton_" + str(i) +
                ".setText(QtCore.QCoreApplication.translate('MainWindow','Delete'))"
            )

            # For updating node - button commands

            exec("self.updateNodeButton_" + str(i) +
                 ".clicked.connect(partial(self.nodeUpdate," + str(i) + "))")

            # For deleting node - button commands

            exec("self.deleteNodeButton_" + str(i) +
                 ".clicked.connect(partial(self.nodeDelete,'" +
                 self.appNodes[i, 0] + "'))")

        # Update plot with new points AFTER listing nodes to prevent user from intercepting plot update
        self.nodePlotUpdate()

    def nodeUpdate(self, rownum):
        #---------------------------------------------------------------------#
        #                       Update node coordinates                       #
        #---------------------------------------------------------------------#
        exec("self._name = self.editNodeName_" + str(rownum) +
             ".text().strip()")
        exec("self._x = self.editXnode_" + str(rownum) + ".text()")
        exec("self._y = self.editYnode_" + str(rownum) + ".text()")
        exec("self._z = self.editZnode_" + str(rownum) + ".text()")
        self._oldname = self.appNodes[rownum, 0]
        self.nodeSubmit(self._name, self._x, self._y, self._z, self._oldname)

    def nodeDelete(self, name):
        #---------------------------------------------------------------------#
        #                             Delete nodes                            #
        #---------------------------------------------------------------------#

        self.nodeExists(name)
        self.truss.RemoveNode(name)
        self.appNodes = delete(self.appNodes, self.check, 0)
        print(f"Deleted Node: {name}")

        # Delete all rows and re-create table
        self.tableClear(self.design.table2editNode)
        self.nodeDisplay()

    def nodePlot(self):
        #---------------------------------------------------------------------#
        #                           Define node plot                          #
        #---------------------------------------------------------------------#
        # Grab necessary variables from design file
        self.widget1plotAreaNode = self.design.widget1plotAreaNode

        # Layout -> Widget
        self.plotAreaNodeLayout = QtWidgets.QVBoxLayout(
            self.widget1plotAreaNode)

        # Plot Canvas
        self.nodeCanvas = MplCanvas(self, dpi=100)

        nodeFig = Figure(dpi=100)
        self.nodeCanvas.node_ax = nodeFig.add_subplot(projection='3d')

        # Plot Toolbar
        self.nodeNavToolbar = NavigationToolbar(self.nodeCanvas,
                                                self.MainWindow)

        # Place canvas and toolbar inside layout
        self.plotAreaNodeLayout.addWidget(self.nodeCanvas)
        self.plotAreaNodeLayout.addWidget(self.nodeNavToolbar)

        # Update which points to plot
        self.nodePlotUpdate()

        # Rename axes
        self.nodeCanvas.node_ax.set_xlabel('X axis')
        self.nodeCanvas.node_ax.set_ylabel('Y axis')
        self.nodeCanvas.node_ax.set_zlabel('Z axis')

        # Show Canvas
        self.nodeCanvas.show()

    def nodePlotUpdate(self):
        #---------------------------------------------------------------------#
        #                          Update Plot Points                         #
        #---------------------------------------------------------------------#

        # Figure settings inside Canvas
        xval = self.appNodes[:, 1].astype(float)
        yval = self.appNodes[:, 2].astype(float)
        zval = self.appNodes[:, 3].astype(float)

        # Clear figure 'node_ax'
        self.nodeCanvas.node_ax.cla()

        # Redefine figure 'node_ax'
        self.nodeCanvas.node_ax.scatter(xval, yval, zval, alpha=0.5)

        # Redraw figure on canvas
        self.nodeCanvas.draw()

    #-------------------------------------------------------------------------#
    #                                                                         #
    #                            Member Functions                             #
    #                                                                         #
    #-------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
#                                                                             #
#                               Initiate Program                              #
#                                                                             #
#-----------------------------------------------------------------------------#

if __name__ == "__main__":
    # Initiate GUI interface
    qapp = QtWidgets.QApplication(argv)
    MainWindow = QtWidgets.QMainWindow()

    # Create 3D Model
    truss = FEModel3D()

    # Initiate program functionality
    app = ApplicationWindow(MainWindow, truss)

    # Show GUI interface
    MainWindow.show()

    # Exit program
    qapp.exec_()
