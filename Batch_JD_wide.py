#!/usr/bin/env python3

#Import Libraries 
import nibabel as nib
import os 
import tkinter as tk
from tkinter import filedialog
from tkinter import font
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.backends.backend_tkagg as tkagg
from scipy.spatial.distance import directed_hausdorff


#Function to convert NIFTI or MGZ into .nii
def mgz_convert(input_mgz,output_nii): 
    try:
        # Load the MGZ file
        img = nib.load(input_mgz)
    
        # Save as NIfTI (.nii)
        nib.save(img, output_nii)
        print("Conversion successful!")
    
    except Exception as e:
        print("Conversion failed:", e)

# Function Converts File format Into nii.seg.nrrd File format 
def nii_seg_convert(input_nii,output_nii_seg):   
    try:
        image = sitk.ReadImage(input_nii)
        sitk.WriteImage(image, output_nii_seg)
        print('Conversion Successful!')

    except Exception as n: 
        print('Conversion Failed:', n)

#Function to auto convert supported file types into nii.seg.nrrd 
def convert_and_generate(input_file1,input_file2):
    supported_extensions = ['.nii.gz', '.mgz', '.nii', '.nii.seg.nrrd']

    if not any(input_file1.lower().endswith(ext) for ext in supported_extensions): 
        print('Error 01: File type not supported for File 1')
        return

    if not any(input_file2.lower().endswith(ext) for ext in supported_extensions):
        print('Error 01: File type not supported for File 2')
        return
    
    subject_names = []
    generated_files = []
    
    if not input_file1.lower().endswith('nii.seg.nrrd'): 
        file_dir_1, filename_1 = os.path.split(input_file1)
        file_root_1, file_ext_1 = os.path.splitext(filename_1)
        file_root_1 = file_root_1[:6]
        subject_names.append(file_root_1)
        
        if input_file1.lower().endswith('nii.gz') or input_file1.lower().endswith('.mgz'):
            output_nii_1 = os.path.join(file_dir_1,file_root_1) + '_pmod.nii'
            output_nii_seg_1 = os.path.join(file_dir_1,file_root_1) + '_pmod.nii.seg.nrrd'
            mgz_convert(input_file1, output_nii_1)
            nii_seg_convert(output_nii_1,output_nii_seg_1)
            generated_files.append(output_nii_seg_1)
    
        elif input_file1.lower().endswith('.nii'):
            output_nii_seg_1 = os.path.join(file_dir_1,file_root_1) + '_pmod.nii.seg.nrrd'
            nii_seg_convert(input_file1, output_nii_seg_1)
            generated_files.append(output_nii_seg_1)
    
    if not input_file2.lower().endswith('.nii.seg.nrrd'):
        file_dir_2, filename_2 = os.path.split(input_file2)
        file_root_2, file_ext_2 = os.path.splitext(filename_2)
        file_root_2 = file_root_2[:6]
        subject_names.append(file_root_2)

        if input_file2.lower().endswith('.nii.gz') or input_file2.lower().endswith('.mgz'):
            output_nii_2 = os.path.join(file_dir_2,file_root_2) + '_smod.nii'
            output_nii_seg_2 = os.path.join(file_dir_2,file_root_2) + '_smod.nii.seg.nrrd'
            mgz_convert(input_file2, output_nii_2)
            nii_seg_convert(output_nii_2, output_nii_seg_2)
            generated_files.append(output_nii_seg_2)

        elif input_file2.lower().endswith('.nii'):
            output_nii_seg_2 = os.path.join(file_dir_2,file_root_2) + '_smod.nii.seg.nrrd'
            nii_seg_convert(input_file2, output_nii_seg_2)
            generated_files.append(output_nii_seg_2)
    
    if output_nii_seg_1 not in generated_files:
        generated_files.append(output_nii_seg_1)

    # Check if the generated files are not already added
    if output_nii_seg_2 not in generated_files:
        generated_files.append(output_nii_seg_2)

    # Ensure unique subject names
    if file_root_1 not in subject_names:
        subject_names.append(file_root_1)

    if file_root_2 not in subject_names:
        subject_names.append(file_root_2)
    return generated_files, subject_names
        

#Function for Segmentation Extraction
def extract_segmentation_info(file_path):
    try: 
        image = sitk.ReadImage(file_path)
        segmentation_array = sitk.GetArrayFromImage(image)
        return segmentation_array
    except Exception as f: 
        print(f'Error 02: Extraction failed for {file_path}: {f}')

#Function to generate binary masks based on user input label
def generate_binary_masks(seg_array, label): 
    try: 
        binary_mask = (seg_array == label).astype(np.uint)
        return binary_mask
    except Exception as g: 
        print(f'Error 03: Mask Generation Failed: {g}')
        return None

def calculate_volume(mask):
    voxel_volume = 0.7  * 0.7 * 0.7 
    volume = np.sum(mask) * voxel_volume

    return volume
    
#Dice Overlap Function 
def dice_coeff(mask1,mask2): 
    intersection = np.logical_and(mask1,mask2)
    dice = 2.0 * np.sum(intersection) / (np.sum(mask1) + np.sum(mask2)) 
    return dice 

#Jaccard Index Function
def jaccard_index(mask1, mask2): 
    intersection = np.logical_and(mask1,mask2)
    union = np.logical_or(mask1,mask2)
    jaccard = np.sum(intersection)/np.sum(union)
    return jaccard

def hausdorff_distance(mask1,mask2): 
    points_mask1 = np.column_stack(np.where(mask1))
    points_mask2 = np.column_stack(np.where(mask2))

    distance1t2 = directed_hausdorff(points_mask1,points_mask2)[0]
    distance2t1 = directed_hausdorff(points_mask2,points_mask1)[0]

    #returns maximum distance
    hausdorff = max(distance1t2,distance2t1)

    return hausdorff

def update_display(subject_name, dice, jaccard, vol_1, vol_2, hausdorff): 
    display_text = f"Subject: {subject_name}\nDice Coeff.: {dice}\nJaccard Index: {jaccard}\nVolume 1: {vol_1}mm^3\nVolume 2: {vol_2}mm^3\nHausdorff Distance: {hausdorff}\n\n"
    current_text = display_panel.get("1.0", tk.END)  # Get current text in the panel
    display_panel.delete("1.0", tk.END)  # Clear the current content of the panel
    display_panel.insert(tk.END, current_text + display_text)  # Append new content
    display_panel.see(tk.END)  # Scroll to the end

# Global variables to store binary masks
binary_mask_1 = None
binary_mask_2 = None
canvas = None
ax1 = None
ax2 = None
# Function to display 3D mask plot
def display_3d_mask_plot():
    global ax1, ax2, canvas

    if binary_mask_1 is not None and binary_mask_2 is not None:
        if canvas is not None:
            canvas.get_tk_widget().destroy()  # Destroy the existing canvas to update with a new one

        fig = plt.figure(figsize=(5, 3))
        ax1 = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122, projection='3d')

        x1, y1, z1 = binary_mask_1.nonzero()
        ax1.scatter(x1, y1, z1, c=z1, cmap= 'magma', marker= ".")
        ax1.axis("off")

        x2, y2, z2 = binary_mask_2.nonzero()
        ax2.scatter(x2, y2, z2, c=z2 , cmap= 'viridis', marker= ".")
        ax2.axis("off")

        ax1.set_xticklabels([])
        ax1.set_yticklabels([])
        ax1.set_zticklabels([])

        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        ax2.set_zticklabels([])
        if canvas is not None and canvas.get_tk_widget().winfo_exists():
            canvas.get_tk_widget().delete('all')
            
        canvas = tkagg.FigureCanvasTkAgg(fig, master=frame_middle)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.grid(row=0, column=0, padx=5, pady=5)
        canvas.draw()
 


#Batch processing function
def batch_processing(directory_path_1, directory_path_2, seg_label_1, seg_label_2):
    global binary_mask_1,binary_mask_2

    files_1 = sorted(os.listdir(directory_path_1))
    files_2 = sorted(os.listdir(directory_path_2))

    if len(files_1) == 0 or len(files_2) == 0:
        print("One or both directories are empty.")
        return
    
    # Skip the first file (.DS_Store in this case) if it exists
    if len(files_1) > 0 and len(files_2) > 0:
        if files_1[0] == '.DS_Store':
            files_1 = files_1[1:]
        if files_2[0] == '.DS_Store':
            files_2 = files_2[1:]

    # Ensure the number of files is the same in both directories 
    num_files = min(len(files_1), len(files_2))

    for i in range(num_files): 
        file_path_1 = os.path.join(directory_path_1, files_1[i])
        file_path_2 = os.path.join(directory_path_2, files_2[i])

        if os.path.isfile(file_path_1) and os.path.isfile(file_path_2): 
            generated_nrrd_files, subject_names = convert_and_generate(file_path_1, file_path_2)
            
            for k, (generated_nrrd_file, subject_name) in enumerate(zip(generated_nrrd_files, subject_names)):
                seg_arrays = []
                for nrrd_file_path in generated_nrrd_files:
                    seg_array = extract_segmentation_info(nrrd_file_path)
                    if seg_array is not None:
                        seg_arrays.append(seg_array)
                    else: 
                        print(f'Error: Segmentation info extraction failed for file: {nrrd_file_path}')
                        continue
                        
                if len(seg_arrays) == 2:
                    seg_array_1 = seg_arrays[0]
                    seg_array_2 = seg_arrays[1]

                    binary_mask_1 = generate_binary_masks(seg_array_1, seg_label_1)
                    binary_mask_2 = generate_binary_masks(seg_array_2, seg_label_2)

                    vol_1 = calculate_volume(binary_mask_1)
                    vol_2 = calculate_volume(binary_mask_2)

                    dice = dice_coeff(binary_mask_1, binary_mask_2)
                    jaccard = jaccard_index(binary_mask_1, binary_mask_2)

                    hausdorff = hausdorff_distance(binary_mask_1, binary_mask_2)

                    if binary_mask_1 is None and binary_mask_2 is None:
                        binary_mask_1 = generate_binary_masks(seg_array_1, seg_label_1)
                        binary_mask_2 = generate_binary_masks(seg_array_2, seg_label_2)
                
                # Display 3D mask plot for the first pair of files
                        

            update_display(subject_name, dice, jaccard, vol_1, vol_2, hausdorff)


def batch_process(): 
    directory_path_1 = entry_path1.get()
    directory_path_2 = entry_path2.get()
    seg_label_1 = int(entry_label1.get())
    seg_label_2 = int(entry_label2.get())

    
    batch_processing(directory_path_1,directory_path_2,seg_label_1,seg_label_2)
    display_3d_mask_plot()

def browse_button(entry): 
    filename = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, filename)


def set_style(widget):
    widget.config(font=("Arial", 10))  # Modify widget's font
    if isinstance(widget, tk.Button):  # Configure button styles
        widget.config(fg="white", bg="#4CAF50")  # Set foreground and background color for buttons

# GUI Setup
root = tk.Tk()
root.title("Batch JD - Roch v1.0")

# Create frames for better organization
frame_left = tk.Frame(root, highlightbackground="black", highlightthickness=1)
frame_left.grid(row=0, column=0, padx=10, pady=10)
frame_left.grid_rowconfigure(0, weight=1)
frame_left.grid_columnconfigure(0, weight=1)

frame_middle = tk.Frame(root, relief=tk.RAISED, borderwidth=2)
frame_middle.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
frame_middle.grid_rowconfigure(0, weight=1)
frame_middle.grid_columnconfigure(0, weight=1)

frame_right = tk.Frame(root)
frame_right.grid(row=0, column=2, padx=10, pady=10)

frame_bottom = tk.Frame(root)
frame_bottom.grid(row=1, column=0, columnspan=3, padx=10, pady=10)
frame_bottom.grid_rowconfigure(0, weight=1)
frame_bottom.grid_columnconfigure(0, weight=1)

# Row 0 - Directory Path 1
label_path1 = tk.Label(frame_bottom, text="Directory Path 1:")
label_path1.grid(row=0, column=0, pady=5)
entry_path1 = tk.Entry(frame_bottom, width=20)
entry_path1.grid(row=0, column=1, padx=5)
button_path1 = tk.Button(frame_bottom, text="Browse", command=lambda: browse_button(entry_path1))
button_path1.grid(row=0, column=2, padx=5)

# Row 1 - Empty space between entries
frame_bottom.rowconfigure(1, weight=1)

# Row 2 - Directory Path 2
label_path2 = tk.Label(frame_bottom, text="Directory Path 2:")
label_path2.grid(row=2, column=0, pady=5)
entry_path2 = tk.Entry(frame_bottom, width=20)
entry_path2.grid(row=2, column=1, padx=5)
button_path2 = tk.Button(frame_bottom, text="Browse", command=lambda: browse_button(entry_path2))
button_path2.grid(row=2, column=2, padx=5)

# Row 3 - Increased spacing between Path 2 and Segmentation Label 1
frame_left.rowconfigure(3, minsize=50)

# Row 4 - Segmentation Label 1
label_label1 = tk.Label(frame_bottom, text="Segmentation Label 1:")
label_label1.grid(row=4, column=0, pady=5)
entry_label1 = tk.Entry(frame_bottom)
entry_label1.grid(row=4, column=1, padx=5)

# Row 5 - Segmentation Label 2
label_label2 = tk.Label(frame_bottom, text="Segmentation Label 2:")
label_label2.grid(row=5, column=0, pady=5)
entry_label2 = tk.Entry(frame_bottom)
entry_label2.grid(row=5, column=1, padx=5)

# Process button
process_button = tk.Button(frame_middle, text="Process", command=batch_process)
process_button.grid(row=2, column=0, padx=5, pady=5, sticky="nsew")
#process_button.configure(bg='black', fg='black')

# Create the canvas for the plot
fig = plt.figure(figsize=(5, 3))
canvas = tkagg.FigureCanvasTkAgg(fig, master=frame_middle)
canvas_widget = canvas.get_tk_widget()
canvas_widget.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")

# Text Widget to display Results
display_panel = tk.Text(frame_middle, height=15, width=30)
display_panel.grid(row=1, column=0, padx=5, pady=5, sticky="nsew")
custom_font = font.Font(family="Courier New", size=10)  # Adjust the size (12 is an example)
display_panel.configure(font=custom_font)

# Run the main GUI loop
root.mainloop()
