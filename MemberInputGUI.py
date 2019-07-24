# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 19:36:17 2018

@author: craig
"""
# Import Tkinter
import tkinter as tk

# Create the root window
MemberGUI = tk.Tk()

# Size the root window
MemberGUI.geometry('293x265')


# Add input boxes
tk.Label(MemberGUI, text = "Member Name:", padx = 5, pady = 5).grid(sticky = "W", row = 0)
tk.Label(MemberGUI, text = "Start Node:", padx = 5, pady = 5).grid(sticky = "W", row = 1)
tk.Label(MemberGUI, text = "End Node:", padx = 5, pady = 5).grid(sticky = "W", row = 2)
tk.Label(MemberGUI, text = "Modulus of Elasticity:", padx = 5, pady = 5).grid(sticky = "W", row = 3)
tk.Label(MemberGUI, text = "Shear Modulus:", padx = 5, pady = 5).grid(sticky = "W", row = 4)
tk.Label(MemberGUI, text = "Moment of Inertia (x-axis):", padx = 5, pady = 5).grid(sticky = "W", row = 5)
tk.Label(MemberGUI, text = "Moment of Inertia (y-axis):", padx = 5, pady = 5).grid(sticky = "W", row = 6)
tk.Label(MemberGUI, text = "Cross-Sectional Area:", padx = 5, pady = 5).grid(sticky = "W", row = 7)
tk.Label(MemberGUI, text = "Torsional Constant:", padx = 5, pady = 5).grid(sticky = "W", row = 8)

tk.Entry(MemberGUI).grid(row = 0, column = 1)
tk.Entry(MemberGUI).grid(row = 1, column = 1)
tk.Entry(MemberGUI).grid(row = 2, column = 1)
tk.Entry(MemberGUI).grid(row = 3, column = 1)
tk.Entry(MemberGUI).grid(row = 4, column = 1)
tk.Entry(MemberGUI).grid(row = 5, column = 1)
tk.Entry(MemberGUI).grid(row = 6, column = 1)
tk.Entry(MemberGUI).grid(row = 7, column = 1)
tk.Entry(MemberGUI).grid(row = 8, column = 1)

# Start the window's main loop
MemberGUI.mainloop()