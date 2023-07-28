#make gif from every image that we have
import cv2
import os                                             
import numpy as np
import imageio

from PIL import Image

image_folder = 'imgs'
video_name = 'result.gif'

images = [img for img in os.listdir(image_folder) if img.endswith(".tga")]
keys = [int(image[5:-4]) for image in images]
xx = np.argsort(keys)
images = np.array(images)[xx]

with imageio.get_writer(video_name, mode='I') as writer:
    for filename in images:
        image = np.array(Image.open(image_folder + "/" + filename))
        writer.append_data(image)
