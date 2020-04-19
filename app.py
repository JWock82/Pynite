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
from numpy import (append, array, delete, empty, float, size, squeeze, vstack,
                   where)
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

        #---------------------------------------------------------------------#
        #                        Application Storage                          #
        #---------------------------------------------------------------------#

        self.appNodes = empty((1, 4))  # A list of the structure's nodes
        self.appMembers = empty((1, 4))  # A list of the structure's members
        self.appMaterials = empty((1, 7))  # A list of materials for reuse
        self.appPlates = empty((1, 5))  # A list of the structure's plates
        self.app__D = {
        }  # A dictionary of the structure's nodal displacements by load combination
        self.appLoadCombos = {
        }  # A dictionary of the structure's load combinations

        self.appNodes = delete(self.appNodes, 0, 0)
        self.appMaterials = delete(self.appMaterials, 0, 0)
        self.appMembers = delete(self.appMembers, 0, 0)

        #---------------------------------------------------------------------#
        #                     Initiate Connect Functions                      #
        #---------------------------------------------------------------------#

        self.nodeStart()
        self.materialStart()
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
        self.design.defineNodeButton.clicked.connect(lambda: self.nodeSubmit(
            self.design.defineNodeName.text(), self.design.defineNodeX.text(),
            self.design.defineNodeY.text(), self.design.defineNodeZ.text()))

        # Initiate plot canvas to allow updating values
        self.nodePlot()

    def materialStart(self):
        # Link member button to 'self.materialSubmit'
        # Inputs are taken after the button is pressed so that the latest values are used.
        self.design.defineMaterialButton.clicked.connect(
            lambda: self.materialSubmit(self.design.defineMaterialName.text(),
                                        self.design.defineMaterialE.text(),
                                        self.design.defineMaterialG.text(),
                                        self.design.defineMaterialIy.text(),
                                        self.design.defineMaterialIz.text(),
                                        self.design.defineMaterialJ.text(),
                                        self.design.defineMaterialA.text()))

    def memberStart(self):
        # Link member button to 'self.memberSubmit'
        # Inputs are taken after the button is pressed so that the latest values are used.
        self.design.defineMemberButton.clicked.connect(
            lambda: self.memberSubmit(
                self.design.defineMemberName.text(),
                self.design.defineMemberINode.currentIndex(),
                self.design.defineMemberJNode.currentIndex(),
                self.design.defineMemberMaterial.currentIndex()))

        # Initiate plot canvas to allow updating values
        self.memberPlot()

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
            widgetToRemove = table.itemAt(i).widget()
            widgetToRemove.setParent(None)
            widgetToRemove.deleteLater()

    def floatEvaluate(self, expr):
        try:
            expr = float(re(eval(expr)))
        except NameError:
            print("Coordinate input must contain valid real numbers.")
            return '', True
        except TypeError:
            print("Coordinate input must contain valid real numbers.")
            return '', True
        except SyntaxError:
            print("Syntax Error.")
            return '', True
        except RecursionError:
            print("Cannot contain recursive input.")
            return '', True
        except (OverflowError, FloatingPointError):
            print("Overflow Error")
            return '', True
        except ZeroDivisionError:
            print("Cannot divide by zero.")
            return '', True
        else:
            return expr, False

    def exists(self, name, array):
        #---------------------------------------------------------------------#
        #                         Check if name exists                        #
        #---------------------------------------------------------------------#

        exists = False
        exec('self._check = squeeze(where(' + array + ' == "' + name +
             '")[0])')
        try:
            int(self._check)
        except TypeError:
            exists = False
        else:
            exists = True
        return exists, self._check

    def duplicate(self, array, find):
        #---------------------------------------------------------------------#
        #                      Check if values are unique                     #
        #---------------------------------------------------------------------#
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

    #-------------------------------------------------------------------------#
    #                                                                         #
    #                             Node Functions                              #
    #                                                                         #
    #-------------------------------------------------------------------------#

    def nodeSubmit(self, name, x, y, z, oldname=""):
        #---------------------------------------------------------------------#
        #                             Define node                             #
        #---------------------------------------------------------------------#
        # Check if input is valid
        if name != "" and x != "" and y != "" and z != "":
            name = name.strip()
            x, errors = self.floatEvaluate(x)
            y, errors = self.floatEvaluate(y)
            z, errors = self.floatEvaluate(z)
            if (errors == True):
                return
        else:
            print("Please enter a value.")
            return

        # The following is for when user edits the displayed text
        if (oldname == ""):
            oldname = name
        else:
            oldname.strip()

        exists, check = self.exists(oldname, 'self.appNodes')

        print("Node name:", name)
        print("x:", x, "m")
        print("y:", y, "m")
        print("z:", z, "m")

        # If exists, delete and replace; else, create new.
        if (exists == True):
            print("Node already exists.")
            self.truss.RemoveNode(oldname)
            self.truss.AddNode(name, x, y, z)

            self.appNodes[check, 0] = name
            self.appNodes[check, 1] = str(x)
            self.appNodes[check, 2] = str(y)
            self.appNodes[check, 3] = str(z)

            # Add members that were deleted
            self.memberUpdateOptions(check)

        if (exists == False):
            # Check if values are used by another element of same type
            find = array([str(x), str(y), str(z)])
            dup = self.duplicate(self.appNodes, find)
            if dup == True:
                print("Nodes must be unique.")
                return

            self.truss.AddNode(name, x, y, z)
            self.appNodes = vstack(
                (self.appNodes, array([name, str(x),
                                       str(y), str(z)])))

        print("-" * 80)

        self.tableClear(self.design.table2updateNode)
        self.nodeDisplay()
        self.nodePlotUpdate()
        self.memberOptions()
        self.tableClear(self.design.table2updateMember)
        self.memberDisplay()
        self.memberPlotUpdate()

    def nodeUpdate(self, rownum):
        #---------------------------------------------------------------------#
        #                             Update node                             #
        #---------------------------------------------------------------------#
        exec("self._name = self.updateNodeName_" + str(rownum) +
             ".text().strip()")
        exec("self._x = self.updateNodeX_" + str(rownum) + ".text()")
        exec("self._y = self.updateNodeY_" + str(rownum) + ".text()")
        exec("self._z = self.updateNodeZ_" + str(rownum) + ".text()")
        self._oldname = self.appNodes[rownum, 0]
        self.nodeSubmit(self._name, self._x, self._y, self._z, self._oldname)

    def nodeDelete(self, name):
        #---------------------------------------------------------------------#
        #                             Delete node                             #
        #---------------------------------------------------------------------#
        exists, check = self.exists(name, 'self.appNodes')
        self.truss.RemoveNode(name)
        self.appNodes = delete(self.appNodes, check, 0)
        print(f"Deleted Node: {name}")

        # Delete all rows and re-create table
        self.tableClear(self.design.table2updateNode)
        self.nodeDisplay()
        self.nodePlotUpdate()
        self.memberDeleteOptions(len(self.appMembers), check, 'node')
        self.tableClear(self.design.table2updateMember)
        self.memberDisplay()
        self.memberOptions()

    def nodeDisplay(self):
        #---------------------------------------------------------------------#
        #                            Display nodes                            #
        #---------------------------------------------------------------------#

        # Label rows
        self.updateNodeSAW = self.design.updateNodeSAW
        self.table2updateNode = self.design.table2updateNode

        self.updateNodeNameLabel = QtWidgets.QLabel(self.updateNodeSAW)
        self.updateNodeNameLabel.setMaximumSize(QtCore.QSize(16777215, 20))
        self.updateNodeNameLabel.setObjectName("updateNodeNameLabel")
        self.table2updateNode.addWidget(self.updateNodeNameLabel, 0, 0, 1, 1)

        self.updateNodeXLabel = QtWidgets.QLabel(self.updateNodeSAW)
        self.updateNodeXLabel.setMaximumSize(QtCore.QSize(16777215, 20))
        self.updateNodeXLabel.setObjectName("updateNodeXLabel")
        self.table2updateNode.addWidget(self.updateNodeXLabel, 0, 1, 1, 1)

        self.updateNodeYLabel = QtWidgets.QLabel(self.updateNodeSAW)
        self.updateNodeYLabel.setMaximumSize(QtCore.QSize(16777215, 20))
        self.updateNodeYLabel.setObjectName("updateNodeYLabel")
        self.table2updateNode.addWidget(self.updateNodeYLabel, 0, 2, 1, 1)

        self.updateNodeZLabel = QtWidgets.QLabel(self.updateNodeSAW)
        self.updateNodeZLabel.setMaximumSize(QtCore.QSize(16777215, 20))
        self.updateNodeZLabel.setObjectName("updateNodeZLabel")
        self.table2updateNode.addWidget(self.updateNodeZLabel, 0, 3, 1, 1)

        # Set label text
        self.updateNodeNameLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "Name"))

        self.updateNodeXLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "Global x (m)"))

        self.updateNodeYLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "Global y (m)"))

        self.updateNodeZLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "Global z (m)"))

        i = 0

        while i <= (len(self.appNodes) - 1):
            row = i + 1
            # node name column
            exec("self.updateNodeName_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateNodeSAW)")
            exec("self.updateNodeName_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(150, 0))")
            exec("self.updateNodeName_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateNodeName_" + str(i) +
                 ".setObjectName('updateNodeName_" + str(i) + "')")
            exec("self.updateNodeName_" + str(i) + ".setText('" +
                 self.appNodes[i, 0] + "')")
            exec("self.table2updateNode.addWidget(self.updateNodeName_" +
                 str(i) + ", " + str(row) + ", 0, 1, 1)")

            # node x column
            exec("self.updateNodeX_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateNodeSAW)")
            exec("self.updateNodeX_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.updateNodeX_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateNodeX_" + str(i) +
                 ".setObjectName('updateNodeX_" + str(i) + "')")
            exec("self.updateNodeX_" + str(i) + ".setText('" +
                 str(self.appNodes[i, 1]) + "')")
            exec("self.table2updateNode.addWidget(self.updateNodeX_" + str(i) +
                 ", " + str(row) + ", 1, 1, 1)")

            # node y column
            exec("self.updateNodeY_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateNodeSAW)")
            exec("self.updateNodeY_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.updateNodeY_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateNodeY_" + str(i) +
                 ".setObjectName('updateNodeY_" + str(i) + "')")
            exec("self.updateNodeY_" + str(i) + ".setText('" +
                 self.appNodes[i, 2] + "')")
            exec("self.table2updateNode.addWidget(self.updateNodeY_" + str(i) +
                 ", " + str(row) + ", 2, 1, 1)")

            # node z column
            exec("self.updateNodeZ_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateNodeSAW)")
            exec("self.updateNodeZ_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.updateNodeZ_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateNodeZ_" + str(i) +
                 ".setObjectName('updateNodeZ_" + str(i) + "')")
            exec("self.updateNodeZ_" + str(i) + ".setText('" +
                 self.appNodes[i, 3] + "')")
            exec("self.table2updateNode.addWidget(self.updateNodeZ_" + str(i) +
                 ", " + str(row) + ", 3, 1, 1)")

            # node update button
            exec("self.updateNodeButton_" + str(i) +
                 " = QtWidgets.QPushButton(self.updateNodeSAW)")
            exec("self.updateNodeButton_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(0, 23))")
            exec("self.updateNodeButton_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(90, 23))")
            exec("self.updateNodeButton_" + str(i) +
                 ".setObjectName('updateNodeButton_" + str(i) + "')")
            exec("self.table2updateNode.addWidget(self.updateNodeButton_" +
                 str(i) + ", " + str(row) + ", 4, 1, 1)")

            # node delete button
            exec("self.deleteNodeButton_" + str(i) +
                 " = QtWidgets.QPushButton(self.updateNodeSAW)")
            exec("self.deleteNodeButton_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(0, 23))")
            exec("self.deleteNodeButton_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(90, 23))")
            exec("self.deleteNodeButton_" + str(i) +
                 ".setObjectName('deleteNodeButton_" + str(i) + "')")
            exec("self.table2updateNode.addWidget(self.deleteNodeButton_" +
                 str(i) + ", " + str(row) + ", 5, 1, 1)")

            # Set button text
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

            i += 1

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

        nodeFig = Figure(dpi=100)
        self.nodeCanvas = FigureCanvasQTAgg(nodeFig)
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
        #                           Update node plot                          #
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
    #                           Materials Functions                           #
    #                                                                         #
    #-------------------------------------------------------------------------#

    def materialSubmit(self, name, E, G, Iy, Iz, J, A, oldname=""):
        # Check if input is valid
        if name != "" and E != "" and G != "" and Iy != "" and Iz != "" and J != "" and A != "":
            name = name.strip()
            E, errors = self.floatEvaluate(E)
            G, errors = self.floatEvaluate(G)
            Iy, errors = self.floatEvaluate(Iy)
            Iz, errors = self.floatEvaluate(Iz)
            J, errors = self.floatEvaluate(J)
            A, errors = self.floatEvaluate(A)
            if (errors == True):
                return
        else:
            print("Please enter a value.")
            return

        # The following is for when user edits the displayed text
        if (oldname == ""):
            oldname = name
        else:
            oldname.strip()

        exists, check = self.exists(oldname, 'self.appMaterials')

        print("Material name:", name)
        print("E:", E, "Pa")
        print("G:", G)
        print("Iy:", Iy)
        print("Iz:", Iz)
        print("J:", J)
        print("A:", A)

        # If exists, delete and replace; else, create new.
        if (exists == True):
            print("Material already exists.")
            self.appMaterials[check, 0] = name
            self.appMaterials[check, 1] = str(E)
            self.appMaterials[check, 2] = str(G)
            self.appMaterials[check, 3] = str(Iy)
            self.appMaterials[check, 4] = str(Iz)
            self.appMaterials[check, 5] = str(J)
            self.appMaterials[check, 6] = str(A)

            # Delete members to update values
            i = 0
            while i <= len(self.appMembers) - 1:
                if self.appMembers[i, 3] == check:
                    self.truss.RemoveMember(self.appMembers[i, 0])
                i += 1
            # Add members that were deleted
            self.memberUpdateOptions(check)
        if (exists == False):
            # Check if values are used by another element of same type
            find = array(
                [str(E),
                 str(G),
                 str(J),
                 str(Iy),
                 str(Iz),
                 str(J),
                 str(A)])
            dup = self.duplicate(self.appMaterials, find)
            if dup == True:
                print("Materials must be unique.")
                return

            self.appMaterials = vstack(
                (self.appMaterials,
                 array(
                     [name,
                      str(E),
                      str(G),
                      str(Iy),
                      str(Iz),
                      str(J),
                      str(A)])))

        print("-" * 80)

        self.tableClear(self.design.table2updateMaterial)
        self.materialDisplay()
        self.memberOptions()
        self.tableClear(self.design.table2updateMember)
        self.memberDisplay()
        self.memberPlotUpdate()

    def materialUpdate(self, rownum):
        #---------------------------------------------------------------------#
        #                     Update material coordinates                     #
        #---------------------------------------------------------------------#
        exec("self._name = self.updateMaterialName_" + str(rownum) +
             ".text().strip()")
        exec("self._E = self.updateMaterialE_" + str(rownum) + ".text()")
        exec("self._G = self.updateMaterialG_" + str(rownum) + ".text()")
        exec("self._Iy = self.updateMaterialIy_" + str(rownum) + ".text()")
        exec("self._Iz = self.updateMaterialIz_" + str(rownum) + ".text()")
        exec("self._J = self.updateMaterialJ_" + str(rownum) + ".text()")
        exec("self._A = self.updateMaterialA_" + str(rownum) + ".text()")
        self._oldname = self.appMaterials[rownum, 0]
        self.materialSubmit(self._name, self._E, self._G, self._Iy, self._Iz,
                            self._J, self._A, self._oldname)

    def materialDelete(self, name):
        #---------------------------------------------------------------------#
        #                           Delete materials                          #
        #---------------------------------------------------------------------#
        exists, check = self.exists(name, 'self.appMaterials')
        self.appMaterials = delete(self.appMaterials, check, 0)
        print(f"Deleted Material: {name}")

        # Delete all rows and re-create table
        self.tableClear(self.design.table2updateMaterial)
        self.materialDisplay()
        self.memberDeleteOptions(len(self.appMembers), check, 'material')
        self.tableClear(self.design.table2updateMember)
        self.memberDisplay()
        self.memberOptions()

    def materialDisplay(self):
        #---------------------------------------------------------------------#
        #                          Display materials                          #
        #---------------------------------------------------------------------#

        # Label rows
        self.updateMaterialSAW = self.design.updateMaterialSAW
        self.table2updateMaterial = self.design.table2updateMaterial

        self.updateMaterialNameLabel = QtWidgets.QLabel(self.updateMaterialSAW)
        self.updateMaterialNameLabel.setObjectName("updateMaterialNameLabel")
        self.table2updateMaterial.addWidget(self.updateMaterialNameLabel, 0, 0,
                                            1, 1)

        self.updateMaterialELabel = QtWidgets.QLabel(self.updateMaterialSAW)
        self.updateMaterialELabel.setObjectName("updateMaterialELabel")
        self.table2updateMaterial.addWidget(self.updateMaterialELabel, 0, 1, 1,
                                            1)

        self.updateMaterialGLabel = QtWidgets.QLabel(self.updateMaterialSAW)
        self.updateMaterialGLabel.setObjectName("updateMaterialGLabel")
        self.table2updateMaterial.addWidget(self.updateMaterialGLabel, 0, 2, 1,
                                            1)

        self.updateMaterialIyLabel = QtWidgets.QLabel(self.updateMaterialSAW)
        self.updateMaterialIyLabel.setObjectName("updateMaterialIyLabel")
        self.table2updateMaterial.addWidget(self.updateMaterialIyLabel, 0, 3,
                                            1, 1)

        self.updateMaterialIzLabel = QtWidgets.QLabel(self.updateMaterialSAW)
        self.updateMaterialIzLabel.setObjectName("updateMaterialIzLabel")
        self.table2updateMaterial.addWidget(self.updateMaterialIzLabel, 0, 4,
                                            1, 1)

        self.updateMaterialJLabel = QtWidgets.QLabel(self.updateMaterialSAW)
        self.updateMaterialJLabel.setObjectName("updateMaterialJLabel")
        self.table2updateMaterial.addWidget(self.updateMaterialJLabel, 0, 5, 1,
                                            1)

        self.updateMaterialALabel = QtWidgets.QLabel(self.updateMaterialSAW)
        self.updateMaterialALabel.setObjectName("updateMaterialALabel")
        self.table2updateMaterial.addWidget(self.updateMaterialALabel, 0, 6, 1,
                                            1)

        # Set label text
        self.updateMaterialNameLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "Name"))

        self.updateMaterialELabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "E"))

        self.updateMaterialGLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "G"))

        self.updateMaterialIyLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "Iy"))

        self.updateMaterialIzLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "Iz"))

        self.updateMaterialJLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "Ix/J"))

        self.updateMaterialALabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "A"))

        i = 0

        while i <= (len(self.appMaterials) - 1):
            row = i + 1
            # material name column
            exec("self.updateMaterialName_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateMaterialSAW)")
            exec("self.updateMaterialName_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(150, 0))")
            exec("self.updateMaterialName_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMaterialName_" + str(i) +
                 ".setObjectName('updateMaterialName_" + str(i) + "')")
            exec("self.updateMaterialName_" + str(i) + ".setText('" +
                 self.appMaterials[i, 0] + "')")
            exec("self.table2updateMaterial.addWidget(self.updateMaterialName_"
                 + str(i) + ", " + str(row) + ", 0, 1, 1)")

            # material E column
            exec("self.updateMaterialE_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateMaterialSAW)")
            exec("self.updateMaterialE_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.updateMaterialE_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMaterialE_" + str(i) +
                 ".setObjectName('updateMaterialE_" + str(i) + "')")
            exec("self.updateMaterialE_" + str(i) + ".setText('" +
                 self.appMaterials[i, 1] + "')")
            exec("self.table2updateMaterial.addWidget(self.updateMaterialE_" +
                 str(i) + ", " + str(row) + ", 1, 1, 1)")

            # material G column
            exec("self.updateMaterialG_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateMaterialSAW)")
            exec("self.updateMaterialG_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.updateMaterialG_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMaterialG_" + str(i) +
                 ".setObjectName('updateMaterialG_" + str(i) + "')")
            exec("self.updateMaterialG_" + str(i) + ".setText('" +
                 self.appMaterials[i, 2] + "')")
            exec("self.table2updateMaterial.addWidget(self.updateMaterialG_" +
                 str(i) + ", " + str(row) + ", 2, 1, 1)")

            # material Iy column
            exec("self.updateMaterialIy_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateMaterialSAW)")
            exec("self.updateMaterialIy_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.updateMaterialIy_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMaterialIy_" + str(i) +
                 ".setObjectName('updateMaterialIy_" + str(i) + "')")
            exec("self.updateMaterialIy_" + str(i) + ".setText('" +
                 self.appMaterials[i, 3] + "')")
            exec("self.table2updateMaterial.addWidget(self.updateMaterialIy_" +
                 str(i) + ", " + str(row) + ", 3, 1, 1)")

            # material Iz column
            exec("self.updateMaterialIz_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateMaterialSAW)")
            exec("self.updateMaterialIz_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.updateMaterialIz_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMaterialIz_" + str(i) +
                 ".setObjectName('updateMaterialIz_" + str(i) + "')")
            exec("self.updateMaterialIz_" + str(i) + ".setText('" +
                 self.appMaterials[i, 4] + "')")
            exec("self.table2updateMaterial.addWidget(self.updateMaterialIz_" +
                 str(i) + ", " + str(row) + ", 4, 1, 1)")

            # material J column
            exec("self.updateMaterialJ_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateMaterialSAW)")
            exec("self.updateMaterialJ_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.updateMaterialJ_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMaterialJ_" + str(i) +
                 ".setObjectName('updateMaterialJ_" + str(i) + "')")
            exec("self.updateMaterialJ_" + str(i) + ".setText('" +
                 self.appMaterials[i, 5] + "')")
            exec("self.table2updateMaterial.addWidget(self.updateMaterialJ_" +
                 str(i) + ", " + str(row) + ", 5, 1, 1)")

            # material A column
            exec("self.updateMaterialA_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateMaterialSAW)")
            exec("self.updateMaterialA_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(10, 0))")
            exec("self.updateMaterialA_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMaterialA_" + str(i) +
                 ".setObjectName('updateMaterialA_" + str(i) + "')")
            exec("self.updateMaterialA_" + str(i) + ".setText('" +
                 self.appMaterials[i, 6] + "')")
            exec("self.table2updateMaterial.addWidget(self.updateMaterialA_" +
                 str(i) + ", " + str(row) + ", 6, 1, 1)")

            # material update button
            exec("self.updateMaterialButton_" + str(i) +
                 " = QtWidgets.QPushButton(self.updateMaterialSAW)")
            exec("self.updateMaterialButton_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(0, 23))")
            exec("self.updateMaterialButton_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(90, 23))")
            exec("self.updateMaterialButton_" + str(i) +
                 ".setObjectName('updateMaterialButton_" + str(i) + "')")
            exec(
                "self.table2updateMaterial.addWidget(self.updateMaterialButton_"
                + str(i) + ", " + str(row) + ", 7, 1, 1)")

            # material delete button
            exec("self.deleteMaterialButton_" + str(i) +
                 " = QtWidgets.QPushButton(self.updateMaterialSAW)")
            exec("self.deleteMaterialButton_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(0, 23))")
            exec("self.deleteMaterialButton_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(90, 23))")
            exec("self.deleteMaterialButton_" + str(i) +
                 ".setObjectName('deleteMaterialButton_" + str(i) + "')")
            exec(
                "self.table2updateMaterial.addWidget(self.deleteMaterialButton_"
                + str(i) + ", " + str(row) + ", 8, 1, 1)")

            # Set button text
            exec(
                "self.updateMaterialButton_" + str(i) +
                ".setText(QtCore.QCoreApplication.translate('MainWindow','Update'))"
            )
            exec(
                "self.deleteMaterialButton_" + str(i) +
                ".setText(QtCore.QCoreApplication.translate('MainWindow','Delete'))"
            )
            # For updating material - button commands
            exec("self.updateMaterialButton_" + str(i) +
                 ".clicked.connect(partial(self.materialUpdate," + str(i) +
                 "))")

            # For deleting material - button commands
            exec("self.deleteMaterialButton_" + str(i) +
                 ".clicked.connect(partial(self.materialDelete,'" +
                 self.appMaterials[i, 0] + "'))")

            i += 1

    #-------------------------------------------------------------------------#
    #                                                                         #
    #                            Member Functions                             #
    #                                                                         #
    #-------------------------------------------------------------------------#

    def memberSubmit(self, name, iNode, jNode, material, oldname=""):
        # Check if input is valid
        if name == "" or iNode == -1 or jNode == -1 or material == -1:
            print("Please enter a value.")
            return
        if iNode == jNode:
            print("Member end-nodes cannot be the same.")
            return

        # The following is for when user edits the displayed text (mistakes)
        if (oldname == ""):
            oldname = name
        else:
            oldname.strip()

        exists, check = self.exists(oldname, 'self.appMembers')

        i = self.appNodes[iNode, 0]
        j = self.appNodes[jNode, 0]
        mat = self.appMaterials[material, 0]

        E = float(self.appMaterials[material, 1])
        G = float(self.appMaterials[material, 2])
        Iy = float(self.appMaterials[material, 3])
        Iz = float(self.appMaterials[material, 4])
        J = float(self.appMaterials[material, 5])
        A = float(self.appMaterials[material, 6])

        print("Member name:", name)
        print("I-node:", i)
        print("J-node:", j)
        print("Material:", mat)

        # If exists, delete and replace; else, create new.
        if (exists == True):
            print("Member already exists.")
            self.truss.RemoveMember(oldname)
            self.truss.AddMember(name, i, j, E, G, Iy, Iz, J, A)
            self.appMembers[check, 0] = name
            self.appMembers[check, 1] = iNode
            self.appMembers[check, 2] = jNode
            self.appMembers[check, 3] = material
        if (exists == False):
            # Check if values are used by another element of same type
            find = array([str(iNode), str(jNode), str(material)])
            dup = self.duplicate(self.appMembers, find)
            if dup == True:
                print("Members must be unique.")
                return

            self.truss.AddMember(name, i, j, E, G, Iy, Iz, J, A)
            self.appMembers = vstack(
                (self.appMembers, array([name, iNode, jNode, material])))

        print("-" * 80)

        self.tableClear(self.design.table2updateMember)

        self.memberDisplay()
        self.memberPlotUpdate()

    def memberUpdateOptions(self, row):
        i = 0
        while i <= len(self.appMembers) - 1:
            name = self.appMembers[i, 0]
            iNode = self.appNodes[int(self.appMembers[i, 1]), 0]
            jNode = self.appNodes[int(self.appMembers[i, 2]), 0]
            E = self.appMaterials[int(self.appMembers[i, 3]), 1].astype(float)
            G = self.appMaterials[int(self.appMembers[i, 3]), 2].astype(float)
            Iy = self.appMaterials[int(self.appMembers[i, 3]), 3].astype(float)
            Iz = self.appMaterials[int(self.appMembers[i, 3]), 4].astype(float)
            J = self.appMaterials[int(self.appMembers[i, 3]), 5].astype(float)
            A = self.appMaterials[int(self.appMembers[i, 3]), 6].astype(float)

            if (self.appMembers[i, 1] == str(row)):
                self.truss.AddMember(name, iNode, jNode, E, G, Iy, Iz, J, A)
            if (self.appMembers[i, 2] == str(row)):
                self.truss.AddMember(name, iNode, jNode, E, G, Iy, Iz, J, A)
            if (self.appMembers[i, 3] == str(row)):
                self.truss.AddMember(name, iNode, jNode, E, G, Iy, Iz, J, A)
            i += 1

    def memberDeleteOptions(self, origlen, row, etype):
        if etype == 'node':
            try:
                i = 0
                while i <= origlen - 1:
                    if (self.appMembers[i, 1] == str(row)) or (
                            self.appMembers[i, 2] == str(row)):
                        self.appMembers = delete(self.appMembers, i, 0)
                        self.memberDeleteOptions(origlen, row, 'node')
                    i += 1
            except IndexError:
                return
        elif etype == 'material':
            try:
                i = 0
                while i <= origlen - 1:
                    if (self.appMembers[i, 3] == str(row)):
                        self.appMembers = delete(self.appMembers, i, 0)
                        self.memberDeleteOptions(origlen, row, 'material')
                    i += 1
            except IndexError:
                return

    def memberOptions(self):
        self.design.defineMemberINode.clear()
        self.design.defineMemberJNode.clear()
        self.design.defineMemberMaterial.clear()

        self.design.defineMemberINode.addItems(self.appNodes[:, 0])
        self.design.defineMemberJNode.addItems(self.appNodes[:, 0])

        self.design.defineMemberMaterial.addItems(self.appMaterials[:, 0])

        self.design.defineMemberINode.setCurrentIndex(-1)
        self.design.defineMemberJNode.setCurrentIndex(-1)
        self.design.defineMemberMaterial.setCurrentIndex(-1)

    def memberUpdate(self, rownum):
        #---------------------------------------------------------------------#
        #                            Update member                            #
        #---------------------------------------------------------------------#
        exec("self._name = self.updateMemberName_" + str(rownum) +
             ".text().strip()")
        exec("self._iNode = self.updateMemberINode_" + str(rownum) +
             ".currentIndex()")
        exec("self._jNode = self.updateMemberJNode_" + str(rownum) +
             ".currentIndex()")
        exec("self._material = self.updateMemberMaterial_" + str(rownum) +
             ".currentIndex()")
        self._oldname = self.appMembers[rownum, 0]
        self.memberSubmit(self._name, self._iNode, self._jNode, self._material,
                          self._oldname)

    def memberDelete(self, name):
        #---------------------------------------------------------------------#
        #                            Delete member                            #
        #---------------------------------------------------------------------#
        exists, check = self.exists(name, 'self.appMembers')
        self.truss.RemoveMember(name)
        self.appMembers = delete(self.appMembers, check, 0)
        print(f"Deleted Member: {name}")

        # Delete all rows and re-create table
        self.tableClear(self.design.table2updateMember)
        self.memberDisplay()
        self.memberPlotUpdate()

    def memberDisplay(self):
        #---------------------------------------------------------------------#
        #                           Display members                           #
        #---------------------------------------------------------------------#

        # Label rows
        self.updateMemberSAW = self.design.updateMemberSAW
        self.table2updateMember = self.design.table2updateMember
        self.updateMemberNameLabel = QtWidgets.QLabel(self.updateMemberSAW)
        self.updateMemberNameLabel.setMaximumSize(QtCore.QSize(16777215, 20))
        self.updateMemberNameLabel.setObjectName("updateMemberNameLabel")
        self.table2updateMember.addWidget(self.updateMemberNameLabel, 0, 0, 1,
                                          1)

        self.updateMemberINodeLabel = QtWidgets.QLabel(self.updateMemberSAW)
        self.updateMemberINodeLabel.setMaximumSize(QtCore.QSize(16777215, 20))
        self.updateMemberINodeLabel.setObjectName("updateMemberINodeLabel")
        self.table2updateMember.addWidget(self.updateMemberINodeLabel, 0, 1, 1,
                                          1)

        self.updateMemberJNodeLabel = QtWidgets.QLabel(self.updateMemberSAW)
        self.updateMemberJNodeLabel.setMaximumSize(QtCore.QSize(16777215, 20))
        self.updateMemberJNodeLabel.setObjectName("updateMemberJNodeLabel")
        self.table2updateMember.addWidget(self.updateMemberJNodeLabel, 0, 2, 1,
                                          1)

        self.updateMemberMaterialLabel = QtWidgets.QLabel(self.updateMemberSAW)
        self.updateMemberMaterialLabel.setMaximumSize(
            QtCore.QSize(16777215, 20))
        self.updateMemberMaterialLabel.setObjectName(
            "updateMemberMaterialLabel")
        self.table2updateMember.addWidget(self.updateMemberMaterialLabel, 0, 3,
                                          1, 1)

        # Set label text
        self.updateMemberNameLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "Name"))
        self.updateMemberINodeLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "I Node"))
        self.updateMemberJNodeLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "J Node"))
        self.updateMemberMaterialLabel.setText(
            QtCore.QCoreApplication.translate("MainWindow", "Material"))

        i = 0

        while i <= (len(self.appMembers) - 1):
            row = i + 1
            # member name column
            exec("self.updateMemberName_" + str(i) +
                 " = QtWidgets.QLineEdit(self.updateMemberSAW)")
            exec("self.updateMemberName_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(150, 0))")
            exec("self.updateMemberName_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMemberName_" + str(i) +
                 ".setObjectName('updateMemberName_" + str(i) + "')")
            exec("self.updateMemberName_" + str(i) + ".setText('" +
                 self.appMembers[i, 0] + "')")
            exec("self.table2updateMember.addWidget(self.updateMemberName_" +
                 str(i) + ", " + str(row) + ", 0, 1, 1)")

            # member i node column
            exec("self.updateMemberINode_" + str(i) +
                 " = QtWidgets.QComboBox(self.updateMemberSAW)")
            exec("self.updateMemberINode_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(20, 0))")
            exec("self.updateMemberINode_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMemberINode_" + str(i) +
                 ".setObjectName('updateMemberINode_" + str(i) + "')")
            exec("self.table2updateMember.addWidget(self.updateMemberINode_" +
                 str(i) + ", " + str(row) + ", 1, 1, 1)")

            # member j node column
            exec("self.updateMemberJNode_" + str(i) +
                 " = QtWidgets.QComboBox(self.updateMemberSAW)")
            exec("self.updateMemberJNode_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(20, 0))")
            exec("self.updateMemberJNode_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMemberJNode_" + str(i) +
                 ".setObjectName('updateMemberJNode_" + str(i) + "')")
            exec("self.table2updateMember.addWidget(self.updateMemberJNode_" +
                 str(i) + ", " + str(row) + ", 2, 1, 1)")

            # member material column
            exec("self.updateMemberMaterial_" + str(i) +
                 " = QtWidgets.QComboBox(self.updateMemberSAW)")
            exec("self.updateMemberMaterial_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(20, 0))")
            exec("self.updateMemberMaterial_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(16777215, 20))")
            exec("self.updateMemberMaterial_" + str(i) +
                 ".setObjectName('updateMemberMaterial_" + str(i) + "')")
            exec("self.table2updateMember.addWidget(self.updateMemberMaterial_"
                 + str(i) + ", " + str(row) + ", 3, 1, 1)")

            # Clear list items
            exec("self.updateMemberINode_" + str(i) + ".clear()")
            exec("self.updateMemberJNode_" + str(i) + ".clear()")

            exec("self.updateMemberMaterial_" + str(i) + ".clear()")

            # Set list items
            exec("self.updateMemberINode_" + str(i) +
                 ".addItems(self.appNodes[:, 0])")
            exec("self.updateMemberJNode_" + str(i) +
                 ".addItems(self.appNodes[:, 0])")

            exec("self.updateMemberMaterial_" + str(i) +
                 ".addItems(self.appMaterials[:, 0])")

            # Set indexes for combobox
            exec("self.updateMemberINode_" + str(i) + ".setCurrentIndex(" +
                 self.appMembers[i, 1] + ")")

            exec("self.updateMemberJNode_" + str(i) + ".setCurrentIndex(" +
                 self.appMembers[i, 2] + ")")

            exec("self.updateMemberMaterial_" + str(i) + ".setCurrentIndex(" +
                 self.appMembers[i, 3] + ")")

            # member update button
            exec("self.updateMemberButton_" + str(i) +
                 " = QtWidgets.QPushButton(self.updateMemberSAW)")
            exec("self.updateMemberButton_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(0, 23))")
            exec("self.updateMemberButton_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(90, 23))")
            exec("self.updateMemberButton_" + str(i) +
                 ".setObjectName('updateMemberButton_" + str(i) + "')")
            exec("self.table2updateMember.addWidget(self.updateMemberButton_" +
                 str(i) + ", " + str(row) + ", 4, 1, 1)")

            # member delete button
            exec("self.deleteMemberButton_" + str(i) +
                 " = QtWidgets.QPushButton(self.updateMemberSAW)")
            exec("self.deleteMemberButton_" + str(i) +
                 ".setMinimumSize(QtCore.QSize(0, 23))")
            exec("self.deleteMemberButton_" + str(i) +
                 ".setMaximumSize(QtCore.QSize(90, 23))")
            exec("self.deleteMemberButton_" + str(i) +
                 ".setObjectName('deleteMemberButton_" + str(i) + "')")
            exec("self.table2updateMember.addWidget(self.deleteMemberButton_" +
                 str(i) + ", " + str(row) + ", 5, 1, 1)")

            exec(
                "self.updateMemberButton_" + str(i) +
                ".setText(QtCore.QCoreApplication.translate('MainWindow','Update'))"
            )
            exec(
                "self.deleteMemberButton_" + str(i) +
                ".setText(QtCore.QCoreApplication.translate('MainWindow','Delete'))"
            )

            # For updating member - button commands

            exec("self.updateMemberButton_" + str(i) +
                 ".clicked.connect(partial(self.memberUpdate," + str(i) + "))")

            # For deleting member - button commands

            exec("self.deleteMemberButton_" + str(i) +
                 ".clicked.connect(partial(self.memberDelete,'" +
                 self.appMembers[i, 0] + "'))")

            i += 1

    def memberPlot(self):
        #---------------------------------------------------------------------#
        #                          Define member plot                         #
        #---------------------------------------------------------------------#
        # Grab necessary variables from design file
        self.widget1plotAreaMember = self.design.widget1plotAreaMember

        # Layout -> Widget
        self.plotAreaMemberLayout = QtWidgets.QVBoxLayout(
            self.widget1plotAreaMember)

        # Plot Canvas

        memberFig = Figure(dpi=100)
        self.memberCanvas = FigureCanvasQTAgg(memberFig)
        self.memberCanvas.member_ax = memberFig.add_subplot(projection='3d')

        # Plot Toolbar
        self.memberNavToolbar = NavigationToolbar(self.memberCanvas,
                                                  self.MainWindow)

        # Place canvas and toolbar inside layout
        self.plotAreaMemberLayout.addWidget(self.memberCanvas)
        self.plotAreaMemberLayout.addWidget(self.memberNavToolbar)

        # Update which points to plot
        self.memberPlotUpdate()

        # Rename axes
        self.memberCanvas.member_ax.set_xlabel('X axis')
        self.memberCanvas.member_ax.set_ylabel('Y axis')
        self.memberCanvas.member_ax.set_zlabel('Z axis')

        # Show Canvas
        self.memberCanvas.show()

    def memberPlotUpdate(self):
        #---------------------------------------------------------------------#
        #                          Update member plot                         #
        #---------------------------------------------------------------------#

        ixval = empty((1))
        iyval = empty((1))
        izval = empty((1))

        jxval = empty((1))
        jyval = empty((1))
        jzval = empty((1))

        ixval = delete(ixval, 0, 0)
        iyval = delete(iyval, 0, 0)
        izval = delete(izval, 0, 0)
        jxval = delete(jxval, 0, 0)
        jyval = delete(jyval, 0, 0)
        jzval = delete(jzval, 0, 0)

        # Figure settings inside Canvas
        for i in range(len(self.appMembers)):
            rownum = int(self.appMembers[i, 1])
            ixval = append(ixval, self.appNodes[rownum, 1]).astype(float)
            iyval = append(iyval, self.appNodes[rownum, 2]).astype(float)
            izval = append(izval, self.appNodes[rownum, 3]).astype(float)

            rownum = int(self.appMembers[i, 2])
            jxval = append(jxval, self.appNodes[rownum, 1]).astype(float)
            jyval = append(jyval, self.appNodes[rownum, 2]).astype(float)
            jzval = append(jzval, self.appNodes[rownum, 3]).astype(float)

        # Clear figure 'member_ax'
        self.memberCanvas.member_ax.cla()

        # Redefine figure 'member_ax'
        self.memberCanvas.member_ax.plot3D((ixval, iyval, izval),
                                           (jxval, jyval, jzval),
                                           alpha=0.5)

        # Redraw figure on canvas
        self.memberCanvas.draw()


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
