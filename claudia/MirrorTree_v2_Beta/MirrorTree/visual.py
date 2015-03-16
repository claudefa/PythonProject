#!/usr/bin/python
# -*- coding: utf-8 -*-

from tkinter import Frame, Tk, BOTH, Text, Menu, END
import tkinter.filedialog 
from tkinter.messagebox import *
from tkinter import *
from MirrorTree.functions import *

class Example(Frame):
  
    def __init__(self, parent):
        Frame.__init__(self, parent)   
         
        self.parent = parent        
        self.initUI()
        
    def initUI(self):
      
        self.parent.title("MirrorTree 2.0")
        self.pack(fill=BOTH, expand=1)
        
        menubar = Menu(self.parent)
        self.parent.config(menu=menubar)
        
        fileMenu = Menu(menubar)

        submenu = Menu(fileMenu)
        submenu.add_command(label="FASTA file", command=self.onOpenFasta)
        submenu.add_command(label="CLUSTAL align", command=self.onOpenClustal)
        fileMenu.add_cascade(label='Open', menu=submenu, underline=0)
        
        fileMenu.add_separator()

        fileMenu.add_command(label="Save", underline=0, command='')        

        fileMenu.add_separator()

        fileMenu.add_command(label="Exit", underline=0, command=self.callback)
        menubar.add_cascade(label="Menu", underline=0, menu=fileMenu)
        menubar.add_separator()              
        menubar.add_command(label="Execute", underline=0,activeforeground="green", command=self.execute)

        self.txt = Text(self)
        self.txt.pack(fill=BOTH, expand=1)


    def onOpenFasta(self):
      
        ftypes = [('FASTA files', '*.fa'), ('FASTA files', '*.fasta'), ('All files', '*')]
        dlg = tkinter.filedialog.Open(self, filetypes = ftypes)
        fl = dlg.show()
        
        if fl != '':
            text = self.readFile(fl)
            self.txt.insert(END, "FASTA file loaded!")
        return text

    def onOpenClustal(self):
      
        ftypes = [('CLUSTAL files', '*.aln'), ('All files', '*')]
        dlg = tkinter.filedialog.Open(self, filetypes = ftypes)
        fl = dlg.show()
        
        if fl != '':
            text = self.readFile(fl)
            self.txt.insert(END, "ClustalW alignment loaded!\n")
            self.txt.insert(END, "%s\n" %(text))
        return text


    def readFile(self, filename):

        f = open(filename, "r")
        text = f.read()
        return filename

    def callback(self):
        if askokcancel('Warning!', 'Do you really want to quit?'):
            raise SystemExit
        else:
            showinfo('', 'Quit has been cancelled!')


    def onExit(self):
        self.quit()

    def execute(self):
        msg = "\nNot implemented yet!\n"
        self.txt.insert(END, msg)

         

def interface():
  
    root = Tk()
    ex = Example(root)
    root.geometry("500x400+300+300")
    explanation = """Welcome to MirrorTree 2.0! :)"""
    w2 = Label(root, 
               justify=LEFT,
               padx = 10, 
               text=explanation).pack(side="bottom")

    root.mainloop()  


if __name__ == '__main__':
    interface()
