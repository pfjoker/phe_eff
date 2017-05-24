# encoding:utf-8
# @Filename=eff_gui
# @Project=phe_eff
# @Date=2017-05-21.21:13
# @Author=sa
# code is far away from bugs with the god animal protecting
"""
I love animals. They taste delicious.
              ┏┓      ┏┓
            ┏┛┻━━━┛┻┓
            ┃      ☃      ┃
            ┃  ┳┛  ┗┳  ┃
            ┃      ┻      ┃
            ┗━┓      ┏━┛
                ┃      ┗━━━┓
                ┃  神兽保佑    ┣┓
                ┃　永无BUG！   ┏┛
                ┗┓┓┏━┳┓┏┛
                  ┃┫┫  ┃┫┫
                  ┗┻┛  ┗┻┛
"""
__metaclass__ = type
from Tkinter import *
from pho_e import pho_e


class Application(Frame, pho_e):
    """My GUI Application"""

    def __init__(self, master):
        # super(Application, self).__init__(master)
        Frame.__init__(self, master)
        pho_e.__init__(self, "pho_e", "canvasname")
        self.grid()
        self.create_widgets()

    def create_widgets(self):
        self.label1 = Label(self, text="plot eff")
        self.label1.grid(row=0, column=0, sticky=W)
        self.button1 = Button(self, text="draw_origin", command=self.draw_origin)
        self.button1.grid(row=1, column=0, sticky=W)
        self.button2 = Button(self, text="draw_origin", command=lambda: self.draw_mass(0, 3, 0, 5))
        self.button2.grid(row=2, column=0, sticky=W)
        self.button3 = Button(self, text="draw_hHitsFit", command=lambda: self.draw_hHitsFit(0, 3, 0.2, 0.3))
        self.button3.grid(row=3, column=0, sticky=W)
        self.button4 = Button(self, text="draw_hHitsDedx", command=lambda: self.draw_hHitsDedx(0, 3, 0.2, 0.3))
        self.button4.grid(row=4, column=0, sticky=W)
        self.button5 = Button(self, text="draw_hHitsDedx_cut", command=lambda: self.draw_HitsDedx_cut_eff(0, 3))
        self.button5.grid(row=5, column=0, sticky=W)
        self.button6 = Button(self, text="draw_hHitsFit_cut", command=lambda: self.draw_HitsFit_cut_eff(0, 3))
        self.button6.grid(row=6, column=0, sticky=W)
        # self.button2 = Button(self, text="order", command=lambda: self.display())
    def display(self):
        pass
if __name__ == '__main__':
    root = Tk()
    root.title('window')
    root.geometry('500x600')
    app = Application(root)
    app.mainloop()
