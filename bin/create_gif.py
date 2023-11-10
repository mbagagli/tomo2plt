#!/usr/bin/env python

import sys
import glob
import imageio.v2 as imageio

anim_file = '%s/tomo.gif' % sys.argv[1]

frame_rate = 3
frame_duration = 1 / frame_rate

with imageio.get_writer(anim_file, mode='I', duration=frame_duration) as writer:
    filenames = glob.glob(
        './%s/*.png' % sys.argv[1])
    filenames = sorted(filenames)
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
    image = imageio.imread(filename)
    writer.append_data(image)

print("DONE")
