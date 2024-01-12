import cv2
import numpy as np

# fourcc = cv2.VideoWriter_fourcc(*'mp4v') 
frameSize = (640, 480)

out = cv2.VideoWriter('output_video.avi',cv2.VideoWriter_fourcc(*'DIVX'), 30, frameSize)
for j in range(0,2500):
     img = cv2.imread("./u0=2/test"+str(j) + '.png')
     out.write(img)
     print(j)

cv2.destroyAllWindows()
out.release()