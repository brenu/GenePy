from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt
import mplcursors
from decimal import Decimal
import numpy as np
from scipy.signal import find_peaks

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import scrolledtext

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)

fileInput = ""
linesNumber = 10


def ctrlEvent(event):
    if(20 == event.state and event.keysym == 'c'):
        return
    else:
        return "break"


def showInformation(selection):
    selection.annotation.set_text(
        f"{selection.artist.get_label()}\nSignal: {selection.target[1]:.2F}\nNucleotide: {selection.target[0]:.2E}")


def getSequenceInfo(data):
    record = SeqIO.read(fileInput, "abi")

    popup = tk.Toplevel(bg="#eff0f1", width=500, height=200)

    sequenceId = tk.Label(popup, bg="#eff0f1", fg="#555555", font=('bold'), anchor="w", justify="left",
                          text="ID: {}".format(record.id), width=50)
    sequenceName = tk.Label(
        popup, bg="#eff0f1", fg="#555555", font=('bold'), anchor="w", justify="left", text="Name: {}".format(record.name), width=50)
    sequenceDescription = tk.Label(
        popup, bg="#eff0f1", fg="#555555", font=('bold'), anchor="w", justify="left", text="Description: {}".format(record.description), width=50)
    sequenceFeatures = tk.Label(
        popup, bg="#eff0f1", fg="#555555", font=('bold'), anchor="w", justify="left", text="Number of features: {}".format(len(record.features)), width=50)
    moleculeType = tk.Label(
        popup, bg="#eff0f1", fg="#555555", font=('bold'), anchor="w", justify="left", text="Molecule type: {}".format(record.annotations["molecule_type"]), width=50)
    machineModel = tk.Label(
        popup, bg="#eff0f1", fg="#555555", font=('bold'), anchor="w", justify="left", text="Machine model: {}".format(record.annotations["machine_model"].decode("utf-8") if record.annotations["machine_model"] else "Unknown"), width=50)
    runStart = tk.Label(
        popup, bg="#eff0f1", fg="#555555", font=('bold'), anchor="w", justify="left", text="Run start: {}".format(record.annotations["run_start"]), width=50)
    runFinish = tk.Label(
        popup, bg="#eff0f1", fg="#555555", font=('bold'), anchor="w", justify="left", text="Run finish: {}".format(record.annotations["run_finish"]), width=50)
    traceScore = tk.Label(
        popup, bg="#eff0f1", fg="#555555", font=('bold'), anchor="w", justify="left", text="Trace score: {}".format(record.annotations["abif_raw"]["TrSc1"] if "TrSc1" in record.annotations["abif_raw"].keys() else "Unknown"), width=50)
    user = tk.Label(
        popup, bg="#eff0f1", fg="#555555", font=('bold'), anchor="w", justify="left", text="User: {}".format(record.annotations["abif_raw"]["User1"].decode("utf-8") if "User1" in record.annotations["abif_raw"].keys() else "Unknown"), width=50)

    sequenceId.grid(column=0, row=0)
    sequenceName.grid(column=0, row=1)
    sequenceDescription.grid(column=0, row=2)
    sequenceFeatures.grid(column=0, row=3)
    moleculeType.grid(column=0, row=4)
    machineModel.grid(column=0, row=5)
    runStart.grid(column=0, row=6)
    runFinish.grid(column=0, row=7)
    traceScore.grid(column=0, row=8)
    user.grid(column=0, row=9)


def getTextSequence(data):
    record = SeqIO.read(fileInput, "abi")

    popup = tk.Toplevel(bg="#eff0f1", width=500, height=200)

    textSequence = scrolledtext.ScrolledText(
        popup, bg="#eff0f1", width=60, height=25, fg="#444444", state="normal")
    # textSequence.insert('1.0', record.seq)
    textSequence.insert('1.0', record.format("fasta"))
    textSequence.bind("<Key>", lambda e: ctrlEvent(e))

    textSequence.grid(column=0, row=0)


def handleSequence(event):
    record = SeqIO.read(fileInput, "abi")
    lines = linesNumber

    # 9 = G      10 = A        11 = T         12 = C
    channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
    trace = defaultdict(list)
    for c in channels:
        if lines > 1:
            trace[c] = np.array_split(record.annotations["abif_raw"][c], lines)
        else:
            trace[c] = record.annotations["abif_raw"][c]

    fig, axs = plt.subplots(lines)

    for ax in axs:
        ax.format_coord = lambda x, y: ""

    fig.text(0.5, 0.02, 'Nucleotides', ha='center')
    fig.text(0.02, 0.5, 'Signal', va='center', rotation='vertical')

    fig.patch.set_facecolor("#eff0f1")

    for i in range(lines):
        if lines > 1:
            axs[i].plot(trace["DATA9"][i], color="#44BBBB", label="G")
            axs[i].plot(trace["DATA10"][i], color="#FF4455", label="A")
            axs[i].plot(trace["DATA11"][i], color="#4BB543", label="T")
            axs[i].plot(trace["DATA12"][i], color="#FFD540", label="C")
        else:
            axs.plot(trace["DATA9"], color="#44BBBB", label="G")
            axs.plot(trace["DATA10"], color="#FF4455", label="A")
            axs.plot(trace["DATA11"], color="#4BB543", label="T")
            axs.plot(trace["DATA12"], color="#FFD540", label="C")

    handles, labels = plt.gca().get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')

    canvas = FigureCanvasTkAgg(fig,
                               master=right)

    cursor = mplcursors.cursor(hover=True)
    cursor.connect('add', showInformation)

    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=5)

    toolbar = NavigationToolbar2Tk(canvas,
                                   right, pack_toolbar=False)

    toolbar.config(background="#eff0f1")
    toolbar._message_label.config(background="#eff0f1")

    toolbar.toolitems = [t for t in NavigationToolbar2Tk.toolitems if
                         t[0] not in ('Pan',)]

    for button in toolbar.winfo_children():
        button.config(background="#eff0f1")

    toolbar.update()
    toolbar.grid(row=2, column=0, padx=10, pady=5)


window = tk.Tk()
window.configure(bg="#eff0f1")
window.winfo_toplevel().title("GenePy")
window.minsize(800, 500)

left = tk.Frame(window, bg="#eff0f1")
right = tk.Frame(window, bg="#eff0f1")
sep = tk.Frame(window, bg="#3aa1a1", width=1, bd=0)

left.grid(row=0, column=0, padx=0, pady=0)
sep.grid(row=0, column=1, padx=0, pady=15, sticky="ns")
right.grid(row=0, column=2, padx=0, pady=0)

greeting = tk.Label(left, text="Choose a file", bg="#eff0f1")
selectFileButton = tk.Button(
    left,
    text="Select",
    width=25,
    height=2,
    bg="#44BBBB",
    fg="#FFFFFF",
    activebackground="#44AAAA",
    activeforeground="#FFFFFF",
    borderwidth=0,
    relief="sunken",
    highlightthickness=0,
    bd=0
)
sequenceButton = tk.Button(
    left,
    text="Plot Sequence",

    width=25,
    height=2,
    bg="#44BBBB",
    fg="#FFFFFF",
    activebackground="#44AAAA",
    activeforeground="#FFFFFF",
    borderwidth=0,
    relief="sunken",
    highlightthickness=0,
    bd=0
)
getTextSequenceButton = tk.Button(
    left,
    text="Text Sequence",

    width=25,
    height=2,
    bg="#44BBBB",
    fg="#FFFFFF",
    activebackground="#44AAAA",
    activeforeground="#FFFFFF",
    borderwidth=0,
    relief="sunken",
    highlightthickness=0,
    bd=0
)
getSequenceInfoButton = tk.Button(
    left,
    text="Sequence Information",

    width=25,
    height=2,
    bg="#44BBBB",
    fg="#FFFFFF",
    activebackground="#44AAAA",
    activeforeground="#FFFFFF",
    borderwidth=0,
    relief="sunken",
    highlightthickness=0,
    bd=0
)


def handleFileSelect(event):
    global fileInput
    auxFileInput = filedialog.askopenfilename()

    if auxFileInput != "":
        fileInput = auxFileInput
        greeting.configure(text=fileInput.split("/")[-1])
        selectFileButton.configure(text="Select Another")
        selectFileButton.grid(pady=2)
        getTextSequenceButton.grid(row=2, column=0, padx=10, pady=2)
        getSequenceInfoButton.grid(row=3, column=0, padx=10, pady=2)
        sequenceButton.grid(row=4, column=0, padx=10, pady=(
            2, window.winfo_screenheight()-850))


selectFileButton.bind("<Button-1>", handleFileSelect)
sequenceButton.bind("<Button-1>", handleSequence)
getTextSequenceButton.bind("<Button-1>", getTextSequence)
getSequenceInfoButton.bind("<Button-1>", getSequenceInfo)

greeting.grid(row=0, column=0, padx=10, pady=(25, 12))
selectFileButton.grid(row=1, column=0, padx=10,
                      pady=(2, window.winfo_screenheight()-686))

window.mainloop()
