import tkinter as tk
import numpy as np
import random


EMPTY = 0
TREE = 1
FIRE = 2
WATER = 3



WIND_DIRECTIONS = {
    "NONE": (0,0),
    "N": (-1,0),
    "S": (1,0),
    "W": (0,-1),
    "E": (0,1),
    "NW": (-1,-1),
    "NE": (-1,1),
    "SW": (1,-1),
    "SE": (1,1)
}



class ForestFire:

    def __init__(self, size, p=0.01, f=0.001):

        self.size = size
        self.p = p
        self.f = f

        self.wind_direction = "NONE"
        self.wind_speed = 0

        self.grid = np.zeros((size,size), dtype=int)

    def set_wind(self, direction, speed):
        self.wind_direction = direction
        self.wind_speed = speed


    def neighbors(self, x, y):

        for dx in [-1,0,1]:
            for dy in [-1,0,1]:

                if dx == 0 and dy == 0:
                    continue

                nx = x + dx
                ny = y + dy

                if 0 <= nx < self.size and 0 <= ny < self.size:
                    yield nx, ny, dx, dy


    def near_water(self, x, y):

        for nx, ny, _, _ in self.neighbors(x,y):
            if self.grid[nx,ny] == WATER:
                return True

        return False


    def wind_modifier(self, dx, dy):

        wdx, wdy = WIND_DIRECTIONS[self.wind_direction]

        if (dx,dy) == (wdx,wdy):
            return self.wind_speed / 20

        elif (dx,dy) == (-wdx,-wdy):
            return -self.wind_speed / 20

        return 0


    def step(self):

        new_grid = self.grid.copy()

        for x in range(self.size):
            for y in range(self.size):

                cell = self.grid[x,y]

                if cell == WATER:
                    continue

                if cell == FIRE:
                    new_grid[x,y] = EMPTY


                elif cell == TREE:

                    burn = False

                    for nx, ny, dx, dy in self.neighbors(x,y):

                        if self.grid[nx,ny] == FIRE:

                            p0 = 0.5
                            p_mod = p0 + self.wind_modifier(dx,dy)

                            if self.near_water(x,y):
                                p_mod *= 0.5

                            p_mod = max(0,min(1,p_mod))

                            if random.random() < p_mod:
                                burn = True
                                break


                    if burn:
                        new_grid[x,y] = FIRE

                    elif random.random() < self.f:
                        new_grid[x,y] = FIRE


                elif cell == EMPTY:

                    if random.random() < self.p:
                        new_grid[x,y] = TREE


        self.grid = new_grid



CELL_SIZE = 8
SIZE = 80


class FireGUI:

    def __init__(self):

        self.model = ForestFire(SIZE)

        self.root = tk.Tk()
        self.root.title("Forest Fire Model")

        self.canvas = tk.Canvas(
            self.root,
            width=SIZE*CELL_SIZE,
            height=SIZE*CELL_SIZE
        )

        self.canvas.grid(row=0,column=0,columnspan=4)

        tk.Label(self.root,text="Wind direction").grid(row=1,column=0)

        self.wind_dir = tk.StringVar()
        self.wind_dir.set("NONE")

        directions = ["NONE","N","S","W","E","NW","NE","SW","SE"]

        tk.OptionMenu(self.root,self.wind_dir,*directions).grid(row=1,column=1)

        tk.Label(self.root,text="Wind speed").grid(row=1,column=2)

        self.wind_speed = tk.Scale(self.root,from_=0,to=10,orient=tk.HORIZONTAL)
        self.wind_speed.grid(row=1,column=3)


        tk.Button(self.root,text="Start",command=self.start).grid(row=2,column=0)
        tk.Button(self.root,text="Step",command=self.step).grid(row=2,column=1)
        tk.Button(self.root,text="Random forest",command=self.random_forest).grid(row=2,column=2)
        tk.Button(self.root,text="Add water",command=self.add_water).grid(row=2,column=3)

        self.running = False

        self.draw()

        self.root.mainloop()


    def random_forest(self):

        for x in range(SIZE):
            for y in range(SIZE):

                if random.random() < 0.6:
                    self.model.grid[x,y] = TREE
                else:
                    self.model.grid[x,y] = EMPTY

        self.draw()


    def add_water(self):

        for _ in range(200):

            x = random.randint(0,SIZE-1)
            y = random.randint(0,SIZE-1)

            self.model.grid[x,y] = WATER

        self.draw()


    def draw(self):

        self.canvas.delete("all")

        for x in range(SIZE):
            for y in range(SIZE):

                cell = self.model.grid[x,y]

                color = "white"

                if cell == TREE:
                    color = "green"

                elif cell == FIRE:
                    color = "red"

                elif cell == WATER:
                    color = "blue"

                self.canvas.create_rectangle(
                    y*CELL_SIZE,
                    x*CELL_SIZE,
                    (y+1)*CELL_SIZE,
                    (x+1)*CELL_SIZE,
                    fill=color,
                    outline=""
                )


    def step(self):

        self.model.set_wind(
            self.wind_dir.get(),
            self.wind_speed.get()
        )

        self.model.step()
        self.draw()


    def loop(self):

        if self.running:
            self.step()
            self.root.after(80,self.loop)


    def start(self):

        self.running = not self.running

        if self.running:
            self.loop()




FireGUI()
