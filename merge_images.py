import sys
import numpy as np
import PIL
from PIL import Image

def merge_images(img_list, out_file="vertical.png"):


    imgs = [ PIL.Image.open(i) for i in img_list ]
    #min_shape = sorted( [(np.sum(i.size), i.size ) for i in imgs])[0][1]
    imgs_comb = np.vstack([np.asarray(i) for i in imgs])
    imgs_comb = PIL.Image.fromarray(imgs_comb)
    imgs_comb.save(out_file)

if __name__ == "__main__":
    
    img_list = sys.argv[1].split()
    merge_images(img_list)
