
x0 = 1e-5

xtalk_coeffs = {

"OTA33.SCI": 
    ( 
        # Row 0:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 1:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 2:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 3:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 4:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 5:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 6:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 7:
             ((1,x0,x0,x0,x0,x0,x0,x0),
              (x0,1,x0,x0,x0,x0,x0,x0),
              (x0,x0,1,x0,x0,x0,x0,x0),
              (x0,x0,x0,1,x0,x0,x0,x0),
              (x0,x0,x0,x0,1,x0,x0,x0),
              (x0,x0,x0,x0,x0,1,x0,x0),
              (x0,x0,x0,x0,x0,x0,1,x0),
              (x0,x0,x0,x0,x0,x0,x0,1),
              ),
     ), # end OTA33
}

xtalk_coeffs['OTA32.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA31.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA44.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA43.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA42.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA24.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA23.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA22.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA55.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA61.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA16.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA00.SCI'] = xtalk_coeffs['OTA33.SCI']


xtalk_matrix = {}

import numpy

def invert_all_xtalk():
    for ota in xtalk_coeffs:
        print ota
        ota_matrices = []
        
        for row_matrix in xtalk_coeffs[ota]:
            #print row_matrix,"\n\n"
            
            inverted = numpy.linalg.inv(row_matrix)
            ota_matrices.append(inverted)
            
        xtalk_matrix[ota] = ota_matrices


