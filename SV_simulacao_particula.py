# Importando bibliotecas
import os
import numpy as np
import random
from pygame.locals import *
import pygame.display

# Criando e definindo classe para a partícula
class Particle:
    
    """Uma classe estabelecendo uma partícula bidimensional."""
    
    def __init__(self, x, y, r, color, ptype):
        
        """Inicializa as propriedades da partícula, tais como: sua posição, velocidade, raio e cor."""

        self.position = np.array((x, y)) # Posição inicial da partícula em x e em y
        self.velocity = np.array([random.uniform(-100, 100), random.uniform(-100,100)]) # Velocidades em x e em y aleatórias entre -1.0 e 1.0
        self.radius = r # Raio da partícula
        self.acceleration = np.array([0,0]) # Aceleração inicial da partícula igual a 0 em x e em y
        #self.acceleration_old = self.acceleration
        self.mass = 1 # Massa da partícula
        self.color = color # Estilo da partícula
        self.type = ptype # tipo da partícula
            
    # Estabelecendo os Métodos para Particle        
    def updateProperties(self,dt):
        
        """Computa a nova aceleração, posição e velocidade, de acordo com um passo de tempo dt. 
        
        Obs: No momento, está sendo considerada aceleração constante.
        """

        #self.update_acceleration(dt)
        self.update_position(dt)
        self.update_velocity(dt)
        
    #def update_acceleration(self,dt):
        
        #"""Atualiza a aceleração da partícula. Como agora ela é constante, não vai mudar.
        #"""
        
        #self.acceleration_old = self.acceleration
        #self.acceleration = self.acceleration
        
    def update_position(self,dt):
        
        """Atualiza a posição da partícula, de acordo com um passo de tempo dt.
        """

        position_1 = self.position #estabele a posição atual
        position_1 += dt*self.velocity #adiciona nela o quanto andou no período estabelecido e com a velocidade determinada
        #position_2 = self.acceleration*0.5*(dt**2.0) #estabelece o quanto andou com a aceleração
        self.position = position_1 #+ position_2 #soma o movimento causado pela aceleração no proporcionado apenas com a velocidade
        
    def update_velocity(self,dt):
        
        """Atualiza a velocidade da partícula, de acordo com um passo de tempo dt e sua aceleração.
        """

        v_1 = self.velocity #estabele a velocidade atual
        #a_sum = self.acceleration + self.acceleration #soma a aceleração
        #v_2 = a_sum*0.5*dt #termo da velocidade da aceleração
        self.velocity = v_1 #+ v_2 #soma dos termos para obter a velocidade final
        
    def reflect_side(self):
        
        """Situação que a partícula colide com a lateral.
        """
        
        #Inversão da componente velocidade em x
        self.velocity[0] = - self.velocity[0]

    def reflect_top(self):
        
        """Situação em que a partícula colide com a parte superior/inferior.
        """
        
        #Invert the y velocity component
        self.velocity[1] = - self.velocity[1]

    def apply_boundary_cond(self):
        
        """Aplica as condições de contorno. Atualiza a posição da partícula na tentativa de simular condições de contorno periódico, com os contornos da caixa sendo (0, xmax) e (0,ymax).
        """
        
        #Para além com comprimento da caixa
        if (self.position[0] >= xmax):
            self.position[0] = self.position[0] - xmax
        #Antes da caixa começar (origem)
        if (self.position[0] < 0):
            self.position[0] = self.position[0] + xmax
        #Para cima da altura da caixa
        if (self.position[1] >= ymax):
            self.position[1] = self.position[1] - ymax
        #Para baixo do fundo da caixa
        if (self.position[1] < 0):
            self.position[1] = self.position[1] + ymax
            
    def apply_impulse(self,Jx,Jy):

        ''' Compute and apply velocity change due to an impulse
           applied to the particle.

           args:
           Jx - scalar x-component of the impulse vector
           Jy - scalar y-component of the impulse vector 
        '''

        self.velocity[0] = self.velocity[0] + (Jx/self.mass)
        self.velocity[1] = self.velocity[1] + (Jy/self.mass)
        
    # For convenience, map the components of the particle's position and
    # velocity vector onto the attributes x, y, vx and vy.
    @property
    def x(self):
        return self.position[0]
    @x.setter
    def x(self, value):
        self.position[0] = value
    @property
    def y(self):
        return self.position[1]
    @y.setter
    def y(self, value):
        self.position[1] = value
    @property
    def vx(self):
        return self.velocity[0]
    @vx.setter
    def vx(self, value):
        self.velocity[0] = value
    @property
    def vy(self):
        return self.velocity[1]
    @vy.setter
    def vy(self, value):
        self.velocity[1] = value

    def overlaps(self, other):
        
        """Confere se há overlap entre as partículas.
        """
        #*descompacta coordenadas
        dist_centers = np.hypot(*(self.position - other.position))
        if dist_centers < self.radius + other.radius:
            return True
        
        else:
            return False
        

    def draw(self, screen):
        
        """Desenha as partículas na tela."""

        pygame.draw.circle(screen, self.color, (int(self.x), int(self.y)), self.radius)