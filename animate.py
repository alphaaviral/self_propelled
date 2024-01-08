import cv2
import numpy as np

# fourcc = cv2.VideoWriter_fourcc(*'mp4v') 
frameSize = (640, 480)

out = cv2.VideoWriter('output_video.avi',cv2.VideoWriter_fourcc(*'DIVX'), 30, frameSize)
for j in range(0,50):
     img = cv2.imread("./out/test"+str(j) + '.png')
     out.write(img)

cv2.destroyAllWindows()
out.release()