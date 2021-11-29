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

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)

fileInput = ""
linesNumber = 10


def showInformation(selection):
    selection.annotation.set_text(
        f"{selection.artist.get_label()}\nSignal: {selection.target[1]:.2F}\nNucleotide: {selection.target[0]:.2E}")


def getTextSequence(data):
    record = SeqIO.read(fileInput, "abi")

    # 9 = G      10 = A        11 = T         12 = C
    channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
    trace = defaultdict(list)
    for c in channels:
        trace[c] = record.annotations["abif_raw"][c]

    sequenceResult = ""

    valuesDict = {
        9: "G",
        10: "A",
        11: "T",
        12: "C"
    }

    peaks = {}

    peaks["DATA9"], _ = find_peaks(
        record.annotations["abif_raw"]["DATA9"], height=0)
    peaks["DATA10"], _ = find_peaks(
        record.annotations["abif_raw"]["DATA10"], height=0)
    peaks["DATA11"], _ = find_peaks(
        record.annotations["abif_raw"]["DATA11"], height=0)
    peaks["DATA12"], _ = find_peaks(
        record.annotations["abif_raw"]["DATA12"], height=0)

    sequenceRange = len(trace["DATA9"])

    for i in range(sequenceRange):
        for j in range(9, 13):
            if trace["DATA{}".format(j)][i] and trace["DATA{}".format(j)][i] in peaks["DATA{}".format(j)]:
                sequenceResult = sequenceResult + valuesDict[j]
                break

    print(sequenceResult)


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
            axs[i].plot(trace["DATA9"][i], color="#4BB", label="G")
            axs[i].plot(trace["DATA10"][i], color="#f45", label="A")
            axs[i].plot(trace["DATA11"][i], color="#4BB543", label="T")
            axs[i].plot(trace["DATA12"][i], color="#FFD540", label="C")
        else:
            axs.plot(trace["DATA9"], color="#4BB", label="G")
            axs.plot(trace["DATA10"], color="#f45", label="A")
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
    bg="#4bb",
    fg="#fff",
    activebackground="#4aa",
    activeforeground="#fff",
    borderwidth=0,
    relief="sunken",
    highlightthickness=0,
    bd=0
)
sequenceButton = tk.Button(
    left,
    text="Sequence",

    width=25,
    height=2,
    bg="#4bb",
    fg="#fff",
    activebackground="#4aa",
    activeforeground="#fff",
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
    bg="#4bb",
    fg="#fff",
    activebackground="#4aa",
    activeforeground="#fff",
    borderwidth=0,
    relief="sunken",
    highlightthickness=0,
    bd=0
)


def handleFileSelect(event):
    global fileInput
    fileInput = filedialog.askopenfilename()

    if fileInput != "":
        greeting.configure(text=fileInput.split("/")[-1])
        selectFileButton.configure(text="Select Another")
        selectFileButton.grid(pady=2)
        getTextSequenceButton.grid(row=2, column=0, padx=10, pady=2)
        sequenceButton.grid(row=3, column=0, padx=10, pady=(
            2, window.winfo_screenheight()-786))


selectFileButton.bind("<Button-1>", handleFileSelect)
sequenceButton.bind("<Button-1>", handleSequence)
getTextSequenceButton.bind("<Button-1>", getTextSequence)

greeting.grid(row=0, column=0, padx=10, pady=(25, 12))
selectFileButton.grid(row=1, column=0, padx=10,
                      pady=(2, window.winfo_screenheight()-686))

# separator = tk.Frame(
#     window, width=5, height=1, bg="#f55")
# separator.grid(column=2, row=0, padx=10, pady=0, sticky="N")

window.mainloop()
