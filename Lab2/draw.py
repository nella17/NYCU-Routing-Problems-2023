import cv2
import numpy as np
import sys

map_path = sys.argv[1]
img_path = sys.argv[2]

with open(map_path, 'r') as f:
    num_layer = int(f.readline())
    data = [
        [
            [int(x) for x in x.split('/')]
            for x in line.split(' ')
        ]
        for line in f.readlines()
    ]

def color(demand: int, capacity: int, via: int):
    if via == -3:
        color = [0, 0, 0, 0]
    elif via == -1 or via == -2:
        alpha = 2 * demand / max(capacity, 1)
        if alpha < 1:
            color = [255 * (1 - alpha), 255 * alpha, 0, 100]
        elif demand <= capacity:
            color = [0, 255 * (2 - alpha), 255 * (alpha - 1), 200]
        else:
            color = [0, 0, 255, 255]
    else:
        alpha = via / num_layer
        color = [255, 255, 255, 255 * alpha]
    return np.array(color, dtype='uint8')

height, width = len(data), len(data[0])

image = np.zeros((height, width, 4), dtype='uint8')

for i in range(height):
    for j in range(width):
        image[i][j] = color(*data[i][j])

cv2.imwrite(img_path, image)
