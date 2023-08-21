#Importando bibliotecas
import pygame
import numpy as np
from pygame.locals import *
import pygame.display
import argparse
import random
from SV_simulacao_particula import Particle
from SV_simulacao_colisao import Event

#Definindo cores
WHITE=(255,255,255)
RED=(255,0,0)
GREEN=(0,255,0)
BLUE=(0,0,255)
BLACK=(0,0,0)
PURPLE = (128, 0, 128)

#Possibilitando argumentos 'command line'
parser = argparse.ArgumentParser(description='An Event Driven Molecular Dynamics (EDMD) Simulator - Sah-Vi Simulation', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#número de partículas
parser.add_argument('-n', action="store", dest="npart", default=50, type=int, help='number of particles to simulate')
#passo-tempo
parser.add_argument('--dt', action="store", dest="dt", type=float, default=10.0, help='max time-step size between graphical updates - higher is faster')
#largura da caixa
parser.add_argument('-x', action="store", dest="xmax", type=int, default=300.0, help='simulation box width [px]')
#altura da caixa
parser.add_argument('-y', action="store", dest="ymax", type=int, default=300.0, help='simulation box height [px]')
options = parser.parse_args()


def setup(options,time,r):

    """Estabelece as partículas.
    
    Return:
        particles: array com as partículas
    """
    
    color = PURPLE
    particles = initParticles(options.npart,r,options.xmax,options.ymax, color)
    return particles


def simulate(options, particles, time):

    ''' Advances the particle ensemble over the time interval dt, or
       to the next collision time, whichever comes first.  If a
       collision is detected within (time,time+dt) then it's carried
       out and the sim time is updated.
       
       Args:
           options:    a valid EDMD options object
           particles:  lista de objetos da classe Particle
           time:       simualtion time at start of jump

       Returns:
           particles: lista atualizada de objetos da classe Particle
           time:      novo tempo da simulação
    '''
    
    # Computa a lista com tempo das próximas colisões, sejam elas contra a parede ou outra partícula
    coll_list = getCollisionList(particles,options.xmax,options.ymax)
    #Define o passo-tempo fornecido inicialmente
    dt = options.dt
    #Se tiver menos que 1 elemento na lista de colisões
    if (len(coll_list) < 1):
        dt_col = dt + 1 # Atualiza dt em uma unidade no passo-tempo
    #Caso haja elementos
    else:
        #Define a primeira colisão como o tempo
        dt_col = (coll_list[0]).t

    # Checa por colisões no tempo atual
    # Se, agora, for antes da primeira colisão detectada
    if (dt < dt_col):
        # Não há colisões no atual passo-tempo
        #Avança partículas
        particles = advanceParticles(particles, dt)
        #Atualiza o tempo com mais um passo tempo
        time = time + dt
    #Caso seja maior ou igual
    else:
        # Houve colisão no passo tempo dado
        # Então, realiza-a
        # Avança partículas no tempo dt_col
        particles = advanceParticles(particles, dt_col)
        # Define o primeiro evento como o primeiro tempo da lista de colisões
        firstEvent = coll_list[0]
        # Performa a colisão
        particles = performCollision(firstEvent,particles)
        #Atualiza o tempo
        time +=  dt_col
        
    return particles, time    

def getWallCollisionTime(Part,xmax,ymax):
    
    ''' Compute the first collision time with between
       particle Part and any wall 

    returns object attributes
    t  -  first collision time
    wall - wall associated with first collision

     locals vars
    t_side # side wall collision time
    t_ud   # top or bottom wall collision time
    w_side # which side wall ('r' or 'l')
    w_ud   # top or bottom wall first? ('u' ,d')
    '''
    
    """Computa o primeiro tempo de colisão entre partícula (Part) e qualquer parede.
    
    Variáveis locais: 
        t_side: side wall collision time
        t_ud: top or bottom wall collision time
        w_side: which side wall ('r' or 'l')
        w_ud: top or bottom wall first? ('u' ,d')
    
    Args:
        Part: objeto Particle válido
        xmax: dimensão x da caixa/tela
        ymax: dimensão y da caixa/tela
    
    Returns: 
        t: primeiro tempo de colisão
        wall: parede associada à primeira colisão"""

    # Paredes laterais
    #Se a velocida da partícula em x for maior que 0
    if (Part.velocity[0] > 0):
        # t = (dS - coord_x - r)/velocidade_x
        t_side = (xmax - Part.position[0] - Part.radius)/Part.velocity[0]
        # parede da direita
        w_side = 'r'
    #Se a velocida da partícula em x for menor que 0
    elif (Part.velocity[0] < 0):
        t_side = (0 - Part.position[0] + Part.radius)/Part.velocity[0]
        # parede da esquerda
        w_side = 'l'
    # Caso contrário
    else:
        # Partícula não está se movendo na direção x
        t_side = None
        w_side = None
        

    # Paredes superior e inferior
    # Se a velocidade da partícula em y for maior que 0
    if (Part.velocity[1] > 0):
        t_ud = (ymax - Part.position[1] - Part.radius)/Part.velocity[1]
        #parede de baixo
        w_ud = 'd'
    # Se a velocidade da partícula em y for menor que 0
    elif (Part.velocity[1] < 0):
        t_ud = (0 - Part.position[1] + Part.radius)/Part.velocity[1]
        #parede de cima
        w_ud = 'u'
    else:
        # Partícula não está se movendo na direção y
        t_ud = None
        w_ud = None
        
    # Se não irá colidir em x nem em y    
    if (t_side == None and t_ud == None):
        # Patícula é estácionária
        t = None
        wall = None
    # Caso dor bater primeiro em alguma lateral
    elif (t_side <= t_ud):
        t = t_side
        wall = w_side
    # Caso for atingir ou em cima ou em baixo primeiro
    else:
        t = t_ud
        wall = w_ud
        
    return type('obj', (object,),{'t': t, 'wall': wall})



def  getCollisionTime(Part1, Part2):

    """Computa o tempo até a colisão entre a Part1 e a Part2.
    
    Args:
    
    Returns:
    
    Retorna time como None se nenhuma solução para o tempo de colisão for encontrado."""

    #diferença entre os vetores velocidades
    deltaVel = Part1.velocity - Part2.velocity
    #diferença entre os vetores posição
    deltaPos = Part1.position - Part2.position
    # distância mínima -> soma dos raios (encostadas)
    minDist = Part1.radius + Part2.radius
    #produto vetorial do delta velocidade
    a = np.dot(deltaVel, deltaVel)
    # duas vezes o produto vetorial da posição com a velocidade
    b = 2.0*np.dot(deltaPos,deltaVel)
    # produto vetorial da pisção menos a distância mínima ao quadrado
    c = np.dot(deltaPos,deltaPos) - minDist*minDist
    discrim = b*b - 4*a*c

    if ((discrim > 0) and (b < 0)):
        t1 = (-b - (discrim**0.5))/(2*a)
        return t1
    
    return None


def getCollisionList(particles,xmax,ymax):
    
    ''' Returns an array of collision Event objects, ordered by their
    time attribute (smallest to largest, Nones at the end)

    args:
    particles - an array of Particle objects '''
    
    # lista de colisão -> return
    coll_list = []

    # loop atravésdo da lista de partículas
    for i in range(len(particles)):
        # Verifica a colisão com alguma pareda da partícula sendo analisada
        wall_collision = getWallCollisionTime(particles[i],xmax,ymax)
        #Computa o evento de colidir com a parede
        firstEvent = Event('w',wall_collision.t,i,None,wall_collision.wall)

        #Para cada partícula seguinte da lista
        for j in range(i+1, len(particles)):
            #Se forem partículas diferentes
            if (i != j):
                col_time = getCollisionTime(particles[i],particles[j])

            # Replace firstEvent if coll time is smaller than current
            # Se for ocorrer primeiro uma colisão com outra partícula, o firstEvent é substituído com tal instante
            # firstEvent.time
                if col_time != None:
                    if (col_time < firstEvent.t):
                        firstEvent = Event('p',col_time,i,j,None)


        # Adiciona à lista de colisão se o evento for válido
        if (firstEvent.t != None):
            coll_list.append(firstEvent)
    

    # Organiza a array de Event, conforme o tempo dos eventos, e a retorna
    coll_list = sorted(coll_list,key=lambda event: event.t)
    
    return coll_list


def performCollision(event,particles):
    
    """Aplica o operador de colisão de acordo com Event.
    
    Args:
        event: um objeto Event válido
        particles: lista das partículas
    
    Returs: lista particles atualizada após aplicar impulso em x e y da(s) partícula(s)
    """
    
    #Caso tenha sido uma colisão com a parede
    if (event.wc_log):
        # Performa a colisão com a parede
        # Se a parede foi a direita ou a esquerda
        if (event.wall == 'r' or event.wall == 'l'):
            particles[event.p1_index].reflect_side()
        # Se a parede foi a de baixo ou a de cima
        elif (event.wall == 'u' or event.wall == 'd'):
            particles[event.p1_index].reflect_top()
        else:
            raise RuntimeError("invalid collison event detected.")
    else:
        # Performa colisão entre duas partículas
        #Calcula o impuso
        J = impulse(particles[event.p1_index], particles[event.p2_index])
        # Aplica o impuso em x e em y de cada uma delas
        particles[event.p1_index].apply_impulse(J[0],J[1])
        particles[event.p2_index].apply_impulse(-J[0],-J[1])
        
    return particles

        
def impulse(Part1,Part2):
    
    """Computa o impulso associado com a colisão partícula-partícula.
    
        J = 2*m1*m2*(dv*dr)/(sigma*(m1+m2))
        
        Args:
            Part1: objeto Particle válido (partícula 1 qualquer)
            Part2: objeto Particle válido (partícula 2 qualquer)
        
        Returns: 
            """
    #diferença das posições pré-colisão das partículas
    dr = Part2.position - Part1.position
    #diferença das velocidades pré-colisão das partículas
    dv = Part2.velocity - Part1.velocity
    #soma dos raios
    sigma = Part1.radius + Part2.radius
    
    hmm = 2*Part1.mass*Part2.mass/(Part1.mass + Part2.mass)
    #dot -> produto vetorial
    J = np.dot(dv,dr)*hmm/sigma
    
    #aplicação do impulso em x e em y
    return [J*dr[0]/sigma,J*dr[1]/sigma]


def advanceParticles(particles,dt):
    
    """Avança o conjunto para frente no tempo em uma trajetória de linha reta.
    
    Args:
        particles: lista com as partículas
        dt: passo-tempo
        
    Returns:
        particles: Lista atualizada com as partículas nas suas novas posições."""
    
    for i in range(len(particles)):
        particles[i].updateProperties(dt)
        
    return particles


def initParticles(n,r,xmax, ymax, color):
    
    """Cria n partículas com raio r e cor color em uma caixa/tela com dimensões (xmax, ymax), de forma que garanta que elas sejam criadas em posições aleatórias dentro dessas dimensões.
    
    Args:
        n: número de partículas a serem criadas
        xmax: dimensão x da caixa/tela
        ymax: dimensão y da caixa/tela
        
    Returns:
        particles: lista com objetos Particle"""
       
    particles = []
            
    x = random.randint(r, xmax - r)
    y = random.randint(r, ymax - r)
    #overlap_init = True
    for np in range(n):
        #while overlap_init == True:
        x = float(random.randint(r, xmax - r))
        y = float(random.randint(r, ymax - r))
        particle = (Particle(x, y , r, color))
        for particle_on_screen in particles:
            if particle_on_screen.overlaps(particle): #== True:
                break

            #elif particle_on_screen.overlaps(particle) == False:
                #particles.append(particle)
                #overlap_init == False
        else: 
            particles.append(particle)
    
    return particles
                
def main(options):

    # define the simulation parameters
    r = 7            # particle radius to use
    time = 10.0        # global simulation time
    dt = 0.1
    paused_log = True # paused indicator bool
    color = PURPLE

    ##################
    ## Screen setup ##
    ##################
    pygame.init()
    #tamanho
    (width, height) = (int(options.xmax), int(options.ymax))
    screen = pygame.display.set_mode((width, height))
    #cor
    background_color = (WHITE)
    #escrita
    pygame.display.set_caption('Sah-Vi Simulation')
    
    pygame.mouse.set_visible(1)
    clock = pygame.time.Clock()

    # set-up the particles    
    particles = setup(options,time,r)

    print("Sah-Vi Simulation initialised with", len(particles), "particles.")
    print("Press ESC or click the 'x' to end the simulation.")
    print("Click anywhere on the screen to pause.")

    #The main run loop
    quit_log = False
    paused_log = False
    while not quit_log:
        clock.tick(60)
        if not paused_log:
            particles, time = simulate(options, particles, time)

        #Handle Input Events
        for event in pygame.event.get():
            if event.type == QUIT:
                quit_log = True
            elif event.type == KEYDOWN and event.key == K_ESCAPE:
                quit_log = True
            elif event.type == MOUSEBUTTONDOWN:
                paused_log = not paused_log
            elif event.type == MOUSEBUTTONUP:
                # fist.unpunch()
                pass
            
        for particle in particles:
            particle.updateProperties(dt)

        screen.fill(background_color)
        
        if not paused_log:
            #Desenha todas as partículas na tela
            for particle in particles:
                particle.draw(screen)
            pygame.display.flip()
            

    pygame.quit()
    
if __name__ == '__main__':
    main(options)
