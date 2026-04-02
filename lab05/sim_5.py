import tkinter as tk
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
import numpy as np


size = [131, 135]
X = [0, 142, 284, 427, 569]
Y = [25, 177, 329]

img = Image.open('Magic_8_baII_answers.jpg')
main_img = Image.open('8_magic_ball.jpg')
main_img.thumbnail((400, 400))


def get_fragment(image, x, y, width, height):
    return image.crop((x, y, x + width, y + height))


# ================= GUI =================

root = tk.Tk()
root.title("Magic 8 Ball")

# конвертация изображений
main_photo = ImageTk.PhotoImage(main_img)

label = tk.Label(root, image=main_photo)
label.pack()


def show_fragment(fragment):
    global label

    photo = ImageTk.PhotoImage(fragment)

    label.config(image=photo)
    label.image = photo  # важно сохранить ссылку


# ================= КНОПКА 1 =================

def yes_no():

    p = np.random.random()

    if p <= 0.5:
        fragment1 = get_fragment(img, X[0], Y[0], size[0], size[1])
    else:
        fragment1 = get_fragment(img, X[1], Y[0], size[0], size[1])

    fragment1 = fragment1.resize((size[0]*3, size[1]*3))
    show_fragment(fragment1)


# ================= КНОПКА 2 =================

def extended():

    p = np.random.random()

    P = [1/15, 2/15, 1/15, 1/15, 1/15, 1/15, 1/15, 1/15, 1/15, 0, 1/15, 1/15, 1/15, 1/15, 1/15]

    sum = 0
    for i in range(15):
        P[i] += sum
        sum = P[i]

    for i in range(15):
        if p <= P[i]:
            a = i + 1
            break

    x = a % 5
    y = a // 5

    fragment2 = get_fragment(img, X[x], Y[y], size[0], size[1])
    fragment2 = fragment2.resize((size[0]*3, size[1]*3))
    show_fragment(fragment2)


# ================= КНОПКИ =================

frame = tk.Frame(root)
frame.pack()

btn1 = tk.Button(frame, text="Да/Нет", command=yes_no, width=20)
btn1.grid(row=0, column=0)

btn2 = tk.Button(frame, text="Расширенный ответ", command=extended, width=20)
btn2.grid(row=0, column=1)


root.mainloop()
