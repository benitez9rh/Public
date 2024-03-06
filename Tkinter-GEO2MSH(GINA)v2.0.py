import os
from tkinter import *               #Import everything from Tkinter
from tkinter import ttk
import tkinter.messagebox           #Import messagebox method
from tkinter import filedialog
import subprocess
'''Create main window, call it root'''
root = Tk()
BROWSE_frame = Frame (root)
BROWSE_frame.grid(row = 1, column = 1, columnspan = 3, sticky = EW, pady = 10, padx = 10)
GEO_frame = Frame(root)
GEO_frame.grid(row = 3, column = 1, pady = 20, padx = 10)
GLI_frame = Frame(root)
GLI_frame.grid(row = 3, column = 3, pady = 20, padx = 10)

''' User Defined Variables'''
epsilon = 5.0e-006

'''Methods'''
#description() function (info box explaining how the app works)
def description():
    tkinter.messagebox.showinfo(title = "Description", message = " Simply use the 'Browse' button to browse to the output .geo file from Gmsh.\n\n Then click 'Format', or alternatively 'File' > 'Format' and a new file will be created in the same directory formatted and ready to input in GINA.")

#browse() function to grab the .geo file and assign respective file path
def browse():
    root.inputfilepath = filedialog.askopenfilename(initialdir="C://Desktop/", title = "Select a file", filetypes = ((".geo files","*.geo"), (".txt files", "*.txt"), ("All files", "*.*") ))
    GEO_Text.delete(0.0, END)
    GEO_Text.insert(0.0, root.inputfilepath)
    try:
        f = open(root.inputfilepath, 'r')
    except:
        tkinter.messagebox.showerror(title = "Error", message = "File not found or path is incorrect")
    finally:
        input_file = open(root.inputfilepath, 'r')
        # read the content of the file line by line 
        data_input = input_file.readlines()
        #Submit input file's contents(.geo) onto the respective text box for review by the user
        GEO_file_Text.delete(0.0, END)
        GEO_file_Text.insert(0.0, data_input)
        
def ginaformat():
    try:
        f = open(root.inputfilepath, 'r')
    except:
        tkinter.messagebox.showerror(title = "Error", message = "File not found or path is incorrect")
    finally:
        #from input file path, reads the path to the folder where the file is located, file extension, file name and assigns respective variables
        inputfilename = root.inputfilepath[root.inputfilepath.rfind("/")+1:]
        path = root.inputfilepath[:root.inputfilepath.rfind("/")]
        input_extension = root.inputfilepath[root.inputfilepath.rfind("."):]
        output_extension = ".gli"
        outputfilepath = root.inputfilepath[:root.inputfilepath.rfind(".")] + output_extension
        
         



        # open file in read mode 
        input_file = open(root.inputfilepath, 'r')
  
        # open other file in write mode 
        output_file  = open(outputfilepath, 'w') 
  
        # read the content of the file line by line 
        data_input = input_file.readlines()

        
        #output_file.write("0\t0\t0\t0\n")
        #For each line in input file that starts with Points, read all numbers except the last (corresponding to characteristic length) which GINA doesn't need, and output to file in a new line separated by tab
        PointFlag = 0
        LineFlag = 0
        for i in data_input:
            if i.startswith("Point"):
                PointFlag += 1
                if PointFlag == 1: #i.e. if it's the first point
                    output_file.write("#POINTS\n")
                if i.rfind("//$NAME") != -1: # if it finds the string
                    pn = i[i.rfind("//")+2:]
                    newstr = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in i[:i.rfind("//")])
                    listOfNumbers = [str(n) for n in newstr.split()]
                    output_file.write( '\t'.join (n for n in listOfNumbers[:len(listOfNumbers)] ) + f"\t{pn}" + "\n" )
                else:
                    newstr = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in i)
                    listOfNumbers = [str(n) for n in newstr.split()]
                    output_file.write( '\t'.join (n for n in listOfNumbers[:len(listOfNumbers)] ) + "\n" ) # Used to be     '\t'.join (n for n in listOfNumbers[:len(listOfNumbers)-1] ) + "\n"   which account for an lc in numeric form. If the file has the word "lc" at the end of each point, the -1 is not necessary or it'll cut the last coordinate.
            if i.startswith("Line("): # The "(" is important because later in the file we have "Line Loop" and if we don't include the parenthesis it will create the Line Loops as Lines, and  we don't want that
                LineFlag += 1
                newstr = ''.join((ch if ch in '0123456789' else ' ') for ch in i) # Because the line number and the points that compose them are all integers, we take out ".-e" from the lookup of the string, unlike we do in points coordinates above
                listOfNumbers = [str(n) for n in newstr.split()][1:] # From [1] because [0] is the Line number
                GEOlinenum = [str(n) for n in newstr.split()][0]
                output_file.write("\n#POLYLINE\n")
                output_file.write("$NAME\n")
                output_file.write(f"L{GEOlinenum}\n")
                output_file.write("$TYPE\n")
                output_file.write("2\n")
                output_file.write("$EPSILON\n")
                output_file.write(f"{epsilon}\n")
                output_file.write("$POINTS\n")
                output_file.write( '\n'.join (n for n in listOfNumbers[:len(listOfNumbers)] ) + "\n" )

                
        output_file.write("\n#STOP\n")
        # output_file.write("Mesh.MshFileVersion = 2.2;")
        #Close output file to submit changes and open again to read
        output_file.close()
        output_file  = open(outputfilepath, 'r')
        data_output = output_file.read()
        #Submit input and output files' contents (.geo and .gli, respectively) onto the respective text boxes for review by the user
        GEO_file_Text.delete(0.0, END)
        GEO_file_Text.insert(0.0, data_input)
        GLI_file_Text.delete(0.0, END)
        GLI_file_Text.insert(0.0, data_output)

        # close all files 
        input_file.close() 
        output_file.close()
        #disable text boxes to edit (As these edits won't be replicated in the files)
        GEO_file_Text.config(state=DISABLED)
        GLI_file_Text.config(state=DISABLED) 

'''Labels'''
#Create labels
GEO_Label = Label(BROWSE_frame, text = ".GEO file", pady = 5)
geo_file_label = Label(GEO_frame, text = ".GEO file content", pady = 5)
gli_file_label = Label(GLI_frame, text = ".GLI file content", pady = 5)

'''Text boxes and scrolls'''
#Create scrollbars
GEOscroll = Scrollbar(GEO_frame)
GLIscroll = Scrollbar(GLI_frame)
#Create text boxes
GEO_Text = Text(BROWSE_frame, width = 60, height=2, state=NORMAL)
GEO_file_Text = Text(GEO_frame, height=20, state=NORMAL, yscrollcommand = GEOscroll.set)
GLI_file_Text = Text(GLI_frame, height=20, state=NORMAL, yscrollcommand = GLIscroll.set)
#Configure scrollbar
GEOscroll.config(command=GEO_file_Text.yview)
GLIscroll.config(command=GLI_file_Text.yview)

'''Buttons'''
BrowseButt= Button(BROWSE_frame, text="Browse .geo file", fg="black", font=("Ariel", 9, "bold"), command=browse)
FormatButt= Button(BROWSE_frame, text="Format", fg="black", font=("Ariel", 9, "bold"), command=ginaformat)

'''Allocate widgets'''
GEO_Label.pack(side = LEFT)
GEO_Text.pack(side = LEFT)
BrowseButt.pack(side = LEFT, fill = Y)
FormatButt.pack(side = BOTTOM, fill = X)

geo_file_label.grid(row = 1, column = 1, sticky = W)
gli_file_label.grid(row = 1, column = 1, sticky = W)
GEO_file_Text.grid(row = 2, column = 1)
GLI_file_Text.grid(row = 2, column = 1)
GEOscroll.grid(row = 1, column = 3, rowspan = 3, sticky = NS)
GLIscroll.grid(row = 1, column = 3, rowspan = 3, sticky = NS)

'''Menu bar'''
menubar = Menu(root)
root.config(menu = menubar)

file_menu = Menu(menubar)
menubar.add_cascade(label = "File", menu = file_menu)
file_menu.add_command(label = "Format", command = ginaformat)
file_menu.add_separator()
file_menu.add_command(label="Exit", command=root.quit)

options_menu = Menu(menubar)
menubar.add_cascade(label = "Options", menu = options_menu)
options_menu.add_command(label = "Description", command = description)

'''Prompt on open'''
#Uncomment below if you wish the description() function (info box explaining how the app works) to be prompted at app opening
#description()

'''Add Window title, geometry and create Window's main loop'''
root.title("Python .Geo file Gmsh to Gina Formater")
root.geometry("1500x500")
root.mainloop()