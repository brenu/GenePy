from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt
import mplcursors
from decimal import Decimal
import numpy as np
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

fileInput = ""


def showInformation(selection):
    selection.annotation.set_text(
        f"{selection.artist.get_label()}\nSignal: {selection.target[1]:.2F}\nNucleotide: {selection.target[0]:.2E}")


def handleSequence(event):
    record = SeqIO.read(fileInput, "abi")
    lines = 10

    # 9 = G      10 = A        11 = T         12 = C
    channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
    trace = defaultdict(list)
    for c in channels:
        if lines > 1:
            trace[c] = np.array_split(record.annotations["abif_raw"][c], lines)
        else:
            trace[c] = record.annotations["abif_raw"][c]

    fig, axs = plt.subplots(lines)

    fig.text(0.5, 0.02, 'Nucleotides', ha='center')
    fig.text(0.02, 0.5, 'Signal', va='center', rotation='vertical')

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

    cursor = mplcursors.cursor(hover=True)
    cursor.connect('add', showInformation)
    plt.show()


window = tk.Tk()
window.minsize(800, 500)

left = tk.Frame(window)
right = tk.Frame(window)
sep = tk.Frame(window, bg="#3aa1a1", width=1, bd=0)

left.grid(row=0, column=0, padx=0, pady=0)
sep.grid(row=0, column=1, padx=0, pady=15, sticky="ns")
right.grid(row=0, column=2, padx=0, pady=0)

greeting = tk.Label(left, text="Choose a file")
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
    relief="sunken"
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
    relief="sunken"
)


def handleFileSelect(event):
    global fileInput
    fileInput = filedialog.askopenfilename()

    if fileInput != "":
        greeting.configure(text=fileInput.split("/")[-1])
        selectFileButton.configure(text="Select Another")
        selectFileButton.grid(pady=2)
        sequenceButton.grid(row=2, column=0, padx=10, pady=(
            2, window.winfo_screenheight()-736))


selectFileButton.bind("<Button-1>", handleFileSelect)
sequenceButton.bind("<Button-1>", handleSequence)

greeting.grid(row=0, column=0, padx=10, pady=(25, 12))
selectFileButton.grid(row=1, column=0, padx=10,
                      pady=(2, window.winfo_screenheight()-686))

# separator = tk.Frame(
#     window, width=5, height=1, bg="#f55")
# separator.grid(column=2, row=0, padx=10, pady=0, sticky="N")

window.mainloop()
