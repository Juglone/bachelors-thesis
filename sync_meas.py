import pyautogui
from time import sleep
from sys import argv

#while True:
#    pos = pyautogui.position()
#    print(pos)

meas_length = int(argv[1])

for i in range(5):
    print(f"{5-i} seconds until first click")
    sleep(1)

pyautogui.click(clicks = 2, interval = meas_length)