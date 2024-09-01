from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtCore import Qt
from OpenGL.GL import *
import OpenGL.GL
from OpenGL import GLU
from OpenGL import GLUT
from OpenGL.arrays import vbo
import numpy as np
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
class genericOpenGLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super(genericOpenGLWidget, self).__init__(parent)
        self.atom_positions = []
        self.camera_pos = np.array([0.0, 0.0, 3.0], dtype=np.float32)
        self.camera_front = np.array([0.0, 0.0, -1.0], dtype=np.float32)
        self.camera_up = np.array([0.0, 1.0, 0.0], dtype=np.float32)
        self.yaw = -90.0
        self.pitch = 0.0
        self.lastX = self.width() / 2
        self.lastY = self.height() / 2
        self.first_mouse = True
        self.camera_speed = 0.05
        self.setMouseTracking(True)
    def animatefun(self, animate):
        xoffset = 1
        yoffset = 1
        i = 1
        while i <10000:
            i += 1

            yaw_rotation = np.array([
                [xoffset+1, 0, xoffset+1],
                [0, 1, 0],
                [-xoffset-1, 0, xoffset+1]
            ])
            pitch_rotation = np.array([
                [1, 0, 0],
                [0, yoffset+1, -yoffset-1],
                [0, yoffset+1, yoffset-1]
            ])
            rotation_matrix = np.dot(pitch_rotation, yaw_rotation)
            if np.any(self.atom_positions):
                ind = -1
                for atom in self.atom_positions:
                    ind += 1
                    new_pos = np.dot(rotation_matrix, self.atom_positions[ind]['position'])
                             
                    self.atom_positions[ind]['position'] = new_pos
        else:
            pass
    def zoomout(self):
        self.camera_pos = np.array([0.0, 0.0, 20.0], dtype=np.float32)
        self.camera_front = np.array([0.0, 0.0, -1.0], dtype=np.float32)
        self.camera_up = np.array([0.0, 1.0, 0.0], dtype=np.float32)
        self.update_camera()

    def zoomin(self):
        self.camera_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        self.camera_front = np.array([0.0, 0.0, -1.0], dtype=np.float32)
        self.camera_up = np.array([0.0, 1.0, 0.0], dtype=np.float32)
        self.update_camera()
        
    def update_camera(self):
        glLoadIdentity()
        view = self.camera_pos + self.camera_front
        gluLookAt(self.camera_pos[0], self.camera_pos[1], self.camera_pos[2],
                  view[0], view[1], view[2],
                  self.camera_up[0], self.camera_up[1], self.camera_up[2])
    def mouseMoveEvent(self, event):
        if self.first_mouse:
            self.lastX = event.x()
            self.lastY = event.y()
            self.first_mouse = False

        xoffset = event.x() - self.lastX
        yoffset = self.lastY - event.y()  
        self.lastX = event.x()
        self.lastY = event.y()
        sensitivity = .01
        xoffset *= sensitivity
        yoffset *= sensitivity

        self.yaw += xoffset
        self.pitch += yoffset
        if event.buttons() & Qt.RightButton:
            z_rotation = np.array([0, 0, yoffset])
            if self.pitch > 89.0:
                self.pitch = 89.0
            if self.pitch < -89.0:
                self.pitch = -89.0
            if np.any(self.atom_positions):
                ind = -1
                for atom in self.atom_positions:
                    ind += 1
                    new_pos =  self.atom_positions[ind]['position'] + z_rotation*self.camera_speed*2
                    self.atom_positions[ind]['position'] = new_pos
        if event.buttons() & Qt.LeftButton:
            yaw_rotation = np.array([
                [np.cos(xoffset), 0, np.sin(xoffset)],
                [0, 1, 0],
                [-np.sin(xoffset), 0, np.cos(xoffset)]
            ])
            pitch_rotation = np.array([
                [1, 0, 0],
                [0, np.cos(yoffset), -np.sin(yoffset)],
                [0, np.sin(yoffset), np.cos(yoffset)]
            ])
            rotation_matrix = np.dot(pitch_rotation, yaw_rotation)
            if np.any(self.atom_positions):
                ind = -1
                for atom in self.atom_positions:
                    ind += 1
                    new_pos = np.dot(rotation_matrix, self.atom_positions[ind]['position'])
                             
                    self.atom_positions[ind]['position'] = new_pos
        #print(self.atom_positions)
        self.update()
    def set_coordinates(self, atom_positions):
        self.atom_positions = atom_positions
        print("attmepting uodates")
        print(atom_positions)
        self.update_camera()
        self.paintGL()
    def initializeGL(self):
        glEnable(GL_DEPTH_TEST)
        glClearDepth(1.0)        # Set the depth buffer to its maximum value
        glDepthFunc(GL_LESS)
        glShadeModel(GL_SMOOTH)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_COLOR_MATERIAL)
        glEnable(GL_BLEND)
        #glDisable(GL_CULL_FACE)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
        glClearColor(0.0, 0.0, 0.0, 0.0)
    def resizeGL(self, w, h):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glViewport(0,0,w,h)
        gluPerspective(50.0, w / h, 0.1, 100.0)
        glMatrixMode(GL_MODELVIEW)

    def paintGL(self):
        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        #glMatrixMode(GL_MODELVIEW)
        #glLoadIdentity()
        #glRotatef(30, 1, 1, 0)
        self.update_camera()
        if np.any(self.atom_positions):
            quadric = gluNewQuadric()
            gluQuadricDrawStyle(quadric, GLU_FILL)
            for atom in self.atom_positions:
                glPushMatrix()
                     # GLU_LINE, GLU_POINT
                letter = atom['element']
                #1-5
                if letter == 'H':
                    glColor3f(1.0,1.0,1.0)
                if letter == 'He':
                    glColor3f(0.85,1.0,1.0)
                if letter == 'Li':
                    glColor3f(0.82,0.5,1.0)
                if letter == 'Be':
                    glColor3f(0.76,1.0,0.0)
                if letter == 'B':
                    glColor3f(1.0,0.71,0.71)
                #6-10
                if letter == 'C':
                    glColor3f(0.55,0.55,0.55)
                if letter == 'N':
                    glColor3f(0.19,0.31,0.97)
                if letter == 'O':
                    glColor3f(1.0,0.05,0.05)
                if letter == 'F':
                    glColor3f(0.56,0.88,0.31)
                if letter == 'Ne':
                    glColor3f(1.0,0.71,0.71) 
                #11-15
                if letter == 'Na':
                    glColor3f(0.67,0.36,0.95)
                if letter == 'Mg':
                    glColor3f(0.54,1.0,0.0)
                if letter == 'Al':
                    glColor3f(0.75,0.63,0.63)
                if letter == 'Si':
                    glColor3f(0.94,0.78,0.63)
                if letter == 'P':
                    glColor3f(1.0,0.50,0.0) 
                #16-20
                if letter == 'S':
                    glColor3f(1.0,1.0,0.19)
                if letter == 'Cl':
                    glColor3f(0.12,0.94,0.12)
                if letter == 'Ar':
                    glColor3f(0.5,0.82,0.89)
                if letter == 'K':
                    glColor3f(0.56,0.25,0.83)
                if letter == 'Ca':
                    glColor3f(0.24,1.0,0.0) 
                #21-25
                if letter == 'Sc':
                    glColor3f(0.90,0.90,0.90)
                if letter == 'Ti':
                    glColor3f(0.75,0.76,0.78)
                if letter == 'V':
                    glColor3f(0.65,0.65,0.67)
                if letter == 'Cr':
                    glColor3f(0.54,0.60,0.78)
                if letter == 'Mn':
                    glColor3f(0.61,0.48,0.78) 
                #26-30
                if letter == 'Fe':
                    glColor3f(0.88,0.4,0.2)
                if letter == 'Co':
                    glColor3f(0.94,0.56,0.63)
                if letter == 'Ni':
                    glColor3f(0.31,0.82,0.31)
                if letter == 'Cu':
                    glColor3f(0.78,0.50,0.2)
                if letter == 'Zn':
                    glColor3f(0.49,0.50,0.69) 
                #31-35
                if letter == 'Ga':
                    glColor3f(0.76,0.56,0.56)
                if letter == 'Ge':
                    glColor3f(0.4, 0.56, 0.56)
                if letter == 'As':
                    glColor3f(0.74, 0.50, 0.89)
                if letter == 'Se':
                    glColor3f(1.0, 0.63, 0.0)
                if letter == 'Br':
                    glColor3f(0.65, 0.16, 0.16) 
                #36-40
                if letter == 'Kr':
                    glColor3f(0.36, 0.72, 0.82)
                if letter == 'Rb':
                    glColor3f(0.44,0.18,0.69)
                if letter == 'Sr':
                    glColor3f(0.0,1.0,0.0)
                if letter == 'Y':
                    glColor3f(0.58,1.0,1.0)
                if letter == 'Zr':
                    glColor3f(0.58,0.88,0.88) 
                #41-45
                if letter == 'Nb':
                    glColor3f(0.45, 0.76, 0.79)
                if letter == 'Mo':
                    glColor3f(0.33, 0.71, 0.71)
                if letter == 'Tc':
                    glColor3f(0.23, 0.62, 0.62)
                if letter == 'Ru':
                    glColor3f(0.14, 0.56, 0.56)
                if letter == 'Rh':
                    glColor3f(0.04, 0.49, 0.55) 
                #46-50
                if letter == 'Pd':
                    glColor3f(0.0, 0.41, 0.52)
                if letter == 'Ag':
                    glColor3f(0.75, 0.75, 0.75)
                if letter == 'Cd':
                    glColor3f(1.0, 0.85, 0.56)
                if letter == 'In':
                    glColor3f(0.65, 0.46, 0.45)
                if letter == 'Sn':
                    glColor3f(0.4, 0.5, 0.5) 
                #51-55
                if letter == 'Sb':
                    glColor3f(0.62, 0.39, 0.71)
                if letter == 'Te':
                    glColor3f(0.83, 0.48, 0.0)
                if letter == 'I':
                    glColor3f(0.58, 0.0, 0.58)
                if letter == 'Xe':
                    glColor3f(0.26, 0.62, 0.69)
                if letter == 'Cs':
                    glColor3f(0.34, 0.09, 0.56) 
                #56-60
                if letter == 'Ba':
                    glColor3f(0.0, 0.79, 0.0)
                if letter == 'La':
                    glColor3f(0.44, 0.83, 1.0)
                if letter == 'Ce':
                    glColor3f(1.0, 1.0, 0.78)
                if letter == 'Pr':
                    glColor3f(0.85, 1.0, 0.78)
                if letter == 'Nd':
                    glColor3f(0.78, 1.0, 0.78) 
                #61-65
                if letter == 'Pm':
                    glColor3f(0.64, 1.0, 0.78)
                if letter == 'Sm':
                    glColor3f(0.56, 1.0, 0.78)
                if letter == 'Eu':
                    glColor3f(0.38, 1.0, 0.78)
                if letter == 'Gd':
                    glColor3f(0.27, 1.0, 0.78)
                if letter == 'Tb':
                    glColor3f(0.19, 1.0, 0.78)
                #66-70
                if letter == 'Dy':
                    glColor3f(0.12, 1.0, 0.78)
                if letter == 'Ho':
                    glColor3f(0.0, 1.0, 0.61)
                if letter == 'Er':
                    glColor3f(0.0, 0.9, 0.46)
                if letter == 'Tm':
                    glColor3f(0.0, 0.83, 0.32)
                if letter == 'Yb':
                    glColor3f(0.0, 0.75, 0.22)
                #71-75
                if letter == 'Lu':
                    glColor3f(0.0, 0.67, 0.14)
                if letter == 'Hf':
                    glColor3f(0.3, 0.76, 1.0)
                if letter == 'Ta':
                    glColor3f(0.3, 0.65, 1.0)
                if letter == 'W':
                    glColor3f(0.13, 0.58, 0.84)
                if letter == 'Re':
                    glColor3f(0.15, 0.49, 0.67)    
                #76-80
                if letter == 'Os':
                    glColor3f(0.15, 0.4, 0.59)
                if letter == 'Ir':
                    glColor3f(0.09, 0.33, 0.53)
                if letter == 'Pt':
                    glColor3f(0.82, 0.82, 0.88)
                if letter == 'Au':
                    glColor3f(1.0, 0.82, 0.14)
                if letter == 'Hg':
                    glColor3f(0.72, 0.72, 0.82)
                #81-85
                if letter == 'Tl':
                    glColor3f(0.65, 0.33, 0.3)
                if letter == 'Pb':
                    glColor3f(0.34, 0.35, 0.38)
                if letter == 'Bi':
                    glColor3f(0.62, 0.31, 0.71)
                if letter == 'Po':
                    glColor3f(0.67, 0.36, 0.0)
                if letter == 'At':
                    glColor3f(0.46, 0.31, 0.27)
                #86-90
                if letter == 'Rn':
                    glColor3f(0.26, 0.51, 0.59)
                if letter == 'Fr':
                    glColor3f(0.26, 0.0, 0.4)
                if letter == 'Ra':
                    glColor3f(0.0, 0.49, 0.0)
                if letter == 'Ac':
                    glColor3f(0.44, 0.67, 0.98)
                if letter == 'Th':
                    glColor3f(0.0, 0.73, 1.0)
                #91-95
                if letter == 'Pa':
                    glColor3f(0.0, 0.63, 1.0)
                if letter == 'U':
                    glColor3f(0.0, 0.56, 1.0)
                if letter == 'Np':
                    glColor3f(0.0, 0.5, 1.0)
                if letter == 'Pu':
                    glColor3f(0.0, 0.42, 1.0)
                if letter == 'Am':
                    glColor3f(0.33, 0.36, 0.95)
                #96-100
                if letter == 'Cm':
                    glColor3f(0.47, 0.36, 0.89)
                if letter == 'Bk':
                    glColor3f(0.54, 0.31, 0.89)
                if letter == 'Cf':
                    glColor3f(0.63, 0.21, 0.83)
                if letter == 'Es':
                    glColor3f(0.7, 0.12, 0.83)
                if letter == 'Fm':
                    glColor3f(0.7, 0.12, 0.73)
                #101-105
                if letter == 'Md':
                    glColor3f(0.7, 0.05, 0.65)
                if letter == 'No':
                    glColor3f(0.74, 0.05, 0.53)
                if letter == 'Lr':
                    glColor3f(0.78, 0.0, 0.4)
                if letter == 'Rf':
                    glColor3f(0.8, 0.0, 0.35)
                if letter == 'Db':
                    glColor3f(0.82, 0.0, 0.31)
                #106-109
                if letter == 'Sg':
                    glColor3f(0.85, 0.0, 0.27)
                if letter == 'Bh':
                    glColor3f(0.88, 0.0, 0.22)
                if letter == 'Hs':
                    glColor3f(0.9, 0.0, 0.18)
                if letter == 'Mt':
                    glColor3f(0.92, 0.0, 0.15)
                glTranslatef(*atom['position'])
                gluSphere(quadric, atom['radius'], 50, 50)
                          
                glPopMatrix()
            print("render updated") 
            gluDeleteQuadric(quadric)
                    
        