import csv
import os

reper = "fichiers_csv"

def export_data_csv(file_name, data):
    if not os.path.exists(reper):
        os.makedirs(reper)
        
    path = os.path.join(reper, file_name)
    with open(path, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["squirmer1.x", "squirmer1.y", 
                         "squirmer2.x", "squirmer2.y", 
                         "squirmer1.orientation", "squirmer2.orientation"
                         "Fl_x_sq1", "Fl_y_sq1", 
                         "Fl_x_sq2", "Fl_y_sq2", 
                         "Dist"])
        for line in data:
            writer.writerow(line)

def read_csv_file(file_name):
    data = []
    path = os.path.join(reper, file_name)
    with open(path, 'r', newline='') as csv_file:
        reader = csv.reader(csv_file)
        #Ignore titles
        next(reader)
        for line in data:
            #Make sure they are floats
            fl_line = [float(val) for val in line]
            data.append(fl_line)

    return data