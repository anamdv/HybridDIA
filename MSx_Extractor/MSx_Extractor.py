from tkinter import *
from tkinter import filedialog
from tkinter.ttk import *
from tkinter.messagebox import showinfo
import sys
import time
import os
from pathlib import Path
from pymsfilereader import MSFileReader

def open_txt():
	global txt_filename
	w.filename=filedialog.askopenfilename(initialdir="/", title="Select a file", filetypes=(("txt files", "*.txt"),))
	txt_filename=w.filename
	my_label=Label(w, text=txt_filename).grid(row=0, column=2)

def open_folder():
	global dir_name
	w.dirname=filedialog.askdirectory()
	dir_name=w.dirname
	my_label=Label(w, text=dir_name).grid(row=1, column=2)

def MSx_extract(dir, file):
	my_label=Label(w, text="Starting analysis...").grid(row=2, column=2)
	n=1
	with open(file) as d:
		pairs = [tuple(map(float, i.split(','))) for i in d]
	mylist=[]
	for filename in os.listdir(dir):
		if filename.endswith(".raw"):
			mylist.append(filename)
		w.update()
	for filename in mylist:
		my_label=Label(w, text="Analysing sample: "+str(n)+"/"+str(len(mylist))).grid(row=2, column=2)
		w.update()
		os.chdir(dir)
		rawfile = MSFileReader(filename)
		ScanNumber=[]
		with open("%s.results.txt" % filename, "w") as f:
				sys.stdout=f
				
				for x in pairs:
					os.chdir(dir)
					ScanNumber=[]
					for i in range(1,rawfile.GetNumSpectra()+1):
						if x==rawfile.GetAllMSOrderData(i)[0][0]:
								ScanNumber.append(i)
								print(x,"\t",i ,"\t",rawfile.GetTrailerExtraForScanNum(i)['Multi Inject Info'])
					scanNumber_to_string=' '.join([str(item) for item in ScanNumber])
					command1='msconvert '+dir+'\\'+filename
					command2=' --filter "scanNumber '+scanNumber_to_string+'"'
					command3=' -o '+dir+'\\'+filename+'mzml_output'
					command4=' --outfile '+filename+''.join(map(str,x))
					os.chdir("C:\\Users\\svg993\\AppData\\Local\\Apps\\ProteoWizard 3.0.21246.16b8e05d 64-bit")
					os.system(command1+command2+command3+command4)
				command5='msconvert '+dir+'\\'+filename+'mzml_output'+"\\*.mzml"
				command6=' --merge'
				os.system(command5+command6+command3)
		n+=1
		f.close()
	w.update()
	my_label=Label(w, text="Done!").grid(row=3, column=2)
	



def doMSxExtract():
	MSx_extract(dir_name, txt_filename)


#if __name__ == "__main__":
import tkinter as tk
w = Tk()
w.title("MSx extractor")
lbl2 = Label(w, text="Upload the IS ENDO file")
lbl2.grid(row=0, column=0)
my_btn=Button(w, text="Open File", command=open_txt).grid(row=0, column=1)
lbl = Label(w, text="Indicate the folder that contains the .raw files")
lbl.grid(row=1, column=0)
my_btn=Button(w, text="Select Dir", command=open_folder).grid(row=1, column=1)
Button(w, text="Go!", command=doMSxExtract).grid(row=2, column=0)
w.mainloop()
