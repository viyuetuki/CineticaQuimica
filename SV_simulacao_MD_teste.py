#Importando bibliotecas
import pygame
import numpy as np
from pygame.locals import *
import pygame.display
import argparse
import random
from SV_simulacao_particula import Particle
from SV_simulacao_colisao import Event
import matplotlib.pyplot as plt
import json
from lmfit import Model

#Definindo cores
WHITE=(255,255,255)
RED=(255,0,0)
GREEN=(0,255,0)
BLUE=(0,0,255)
BLACK=(0,0,0)
PURPLE = (128, 0, 128)
ORANGE = (255, 165, 0)

#Possibilitando argumentos 'command line'
parser = argparse.ArgumentParser(description='An Event Driven Molecular Dynamics (EDMD) Simulator - Sah-Vi Simulation', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#número de partículas
parser.add_argument('-nA', action="store", dest="npartA", default=100, type=int, help='number of particle type A to simulate')
parser.add_argument('-nB', action="store", dest="npartB", default=100, type=int, help='number of particle type B to simulate')
#passo-tempo
parser.add_argument('--dt', action="store", dest="dt", type=float, default=10.0, help='max time-step size between graphical updates - higher is faster')
#largura da caixa
parser.add_argument('-x', action="store", dest="xmax", type=int, default=500.0, help='simulation box width [px]')
#altura da caixa
parser.add_argument('-y', action="store", dest="ymax", type=int, default=500.0, help='simulation box height [px]')
#tipo da partícula
parser.add_argument('-ptype', action="store", dest="ptype", type=str, default='A', help='type of the particle being created')
options = parser.parse_args()

def escrever_txt(lista):
    with open('simulacao_con_vel.txt', 'w', encoding='utf-8') as f:
        counter = 0
        for listvel in lista:
            counter +=1
            f.write(f"{counter}:" + '\n')
            for vel in listvel:
                f.write(str(vel) + ",")

def carregar_txt():
    with open('simulacao_con_vel.txt', 'r', encoding='utf-8') as f:
        return f.readlines()
    
    
def overlap_all(list_particles, xmax, ymax, r, color, ptype):
    
    x = random.uniform(r, xmax - r)
    y = random.uniform(r, ymax - r)
    validation = "ok!"
    particle = (Particle(x, y , r, color, ptype))

    for particle_on_screen in list_particles:
        if particle_on_screen.overlaps(particle):
            validation = "overlaped"
            break   
    
    return particle, validation


def initParticles(nA, nB,r1, r2,xmax, ymax, color1, color2, ptype1, ptype2):
    
    """Cria n partículas com raio r e cor color, do tipo ptype, em uma caixa/tela com dimensões (xmax, ymax), de forma que garanta que elas sejam criadas em posições aleatórias dentro dessas dimensões.
    
    Args:
        n: número de partículas a serem criadas
        r: raio das partículas
        xmax: dimensão x da caixa/tela
        ymax: dimensão y da caixa/tela
        color: cor das partículas
        ptype: tipo da partícula
        
    Returns:
        particles: lista com objetos Particle"""
       
    particles = []
            
    #x = random.randint(r, xmax - r)
    #y = random.randint(r, ymax - r)
    #overlap_init = True
    for npA in range(nA):
        #while overlap_init == True:
        #x1 = random.uniform(r1, xmax - r1)
        #y1 = random.uniform(r1, ymax - r1)
        particle1, validation1 = overlap_all(particles, xmax, ymax, r1, color1, ptype1)
        while validation1 == "overlaped":
            particle1, validation1 = overlap_all(particles, xmax, ymax, r1, color1, ptype1)
        
        particles.append(particle1)
        
        #particle1 = (Particle(x1, y1 , r1, color1, ptype1))
        #for particle_on_screen in particles:
        #    if particle_on_screen.overlaps(particle1): #== True:
        #        break

            #elif particle_on_screen.overlaps(particle) == False:
                #particles.append(particle)
                #overlap_init == False
        #else: 
        #    particles.append(particle1)
            
        #x2 = random.uniform(r2, xmax - r2)
        #y2 = random.uniform(r2, ymax - r2)
    for npB in range(nB):

        particle2, validation2 = overlap_all(particles, xmax, ymax, r2, color2, ptype2)
        while validation1 == "overlaped":
            particle2, validation2 = overlap_all(particles, xmax, ymax, r2, color2, ptype2)
        
        particles.append(particle2)
        #particle2 = (Particle(x2, y2, r2, color2, ptype2))
        #for particle_on_screen in particles:
        #    if particle_on_screen.overlaps(particle2): #== True:
        #        break

        #else: 
        #    particles.append(particle2)
    
    return particles


def setup(options, r1, r2, color1, color2, ptype1, ptype2):

    """Estabelece as partículas.
    
    Return:
        particles: array com as partículas
    """
    
    particles = initParticles(options.npartA,options.npartB,r1, r2,options.xmax,options.ymax, color1, color2, ptype1, ptype2)
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

    velocities_dt = []
    for i in range(len(particles)):
        velocity_dt = np.sqrt((particles[i].vx)**2 + (particles[i].vy)**2)
        velocities_dt.append(velocity_dt)
        particles[i].updateProperties(dt)
        
    return particles, velocities_dt


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
        # Calcula o impuso
        J = impulse(particles[event.p1_index], particles[event.p2_index])
        # Aplica o impuso em x e em y de cada uma delas
        particles[event.p1_index].apply_impulse(J[0],J[1])
        particles[event.p2_index].apply_impulse(-J[0],-J[1])

        P1_ptype = particles[event.p1_index].type
        P2_ptype = particles[event.p2_index].type

        if (P1_ptype == "A" and P2_ptype == "B") or (P2_ptype == "A" and P1_ptype == "B"):
            for particle_in_collision in [particles[event.p1_index], particles[event.p2_index]]:
                if particle_in_collision.type == "A":
                    particle_in_collision.type = "C"
                    particle_in_collision.radius = 3
                    particle_in_collision.color = BLUE

                elif particle_in_collision.type == "B":
                    particle_in_collision.type = "D"
                    particle_in_collision.radius = 8
                    particle_in_collision.color = GREEN
     
    return particles


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
        particles, velocities_dt = advanceParticles(particles, dt)
        #plt.hist(velocities_dt, bins=15)
        #plt.title(f"time {time}")
        #Atualiza o tempo com mais um passo tempo
        time = time + dt
    #Caso seja maior ou igual
    else:
        # Houve colisão no passo tempo dado
        # Então, realiza-a
        # Avança partículas no tempo dt_col
        particles, velocities_dt = advanceParticles(particles, dt_col)
        #plt.hist(velocities_dt, bins=15)
        #plt.title(f"time {time}")
        # Define o primeiro evento como o primeiro tempo da lista de colisões
        firstEvent = coll_list[0]
        # Performa a colisão
        particles = performCollision(firstEvent,particles)
        #Atualiza o tempo
        time +=  dt_col
    

    return particles, time, velocities_dt   


def first_order(x, C, k):
    
    y = C*np.exp(-k*x)

    return y
    


def squared(x, C, k):
    
    #y=(-1/(k*x + (1/C)))
    #y = x*t + 1/C
    #y = C*np.exp(-k*x)
    #y = C/(C*k*x + 1)
    y = 1/(k*x+(1/C))

    return y

def squared_2(x, CA, CB, B, k):
    
    #y=(-1/(k*x + (1/C)))
    #y = x*t + 1/C
    #y = C*np.exp(-k*x)
    #y = C/(C*k*x + 1)
    #y = 1/(k*x+(1/C))
    y = np.exp(k*(CA - CB)*x) * CA * (-CB) * B

    return y
    

def fractional(x, C, k):
    
    y = ((x*k)/2 + C**(1/2))**(2)

    return y
    

def derivada_xy(x_vals, y_vals):
    avg_x_vals = []
    deriv_vals = []
    for i in range(len(x_vals) - 1):
        deriv = ( y_vals[i + 1] - y_vals[i] ) / ( x_vals[i + 1] - x_vals[i] )
        avg_x = ( x_vals[i + 1] + x_vals[i] ) / 2

        avg_x_vals.append(avg_x)
        deriv_vals.append(deriv)

    return avg_x_vals, deriv_vals    

                
def main(options):

    # define the simulation parameters
    r1 = 7            # particle radius to use
    r2 = 10
    time = 10.0        # global simulation time
    dt = 0.1
    paused_log = True # paused indicator bool
    color1 = PURPLE
    color2 = ORANGE
    ptype1 = "A"
    ptype2 = "B"
    velocities2_dt = []
    velocities2_dt_total = []
    #counter_A = 0
    particles_A = []
    #counter_B = 0
    particles_B = []
    #counter_C = 0
    particles_C = []
    #counter_D = 0
    particles_D = []
    instant = 0
    time_pass = []
    file_counter = 0 

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
    particles = setup(options,r1, r2, color1, color2, ptype1, ptype2)

    print("Sah-Vi Simulation initialised with", len(particles), "particles.")
    print("Press ESC or click the 'x' to end the simulation.")
    print("Click anywhere on the screen to pause.")

    #The main run loop
    quit_log = False
    paused_log = False
    while not quit_log:
        clock.tick(60)
        if not paused_log:
            particles, time, velocities_dt = simulate(options, particles, time)
            #plt.hist(velocities_dt, bins=10) (certo!!)
            #plt.show()
            #plt.savefig(f"Histograma - dt{time}.png", dpi = 300)

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

        counter_A = 0
        counter_B = 0
        counter_C = 0
        counter_D = 0

        file_counter = 0


        for particle in particles:
            particle.updateProperties(dt)
            velocity2_dt = np.sqrt((particle.vx)**2 + (particle.vy)**2)
            velocities2_dt.append(velocity2_dt)

            if particle.type == "A":
                counter_A += 1

            elif particle.type == "B":
                counter_B += 1
            
            elif particle.type == "C":
                counter_C += 1

            elif particle.type == "D":
                counter_D += 1
            
            #plt.hist(velocities, bins=50)
            #plt.savefig(f"Histograma - dt{time}.png", dpi = 300)

        instant +=1
        particles_A.append(counter_A)
        particles_B.append(counter_B)
        particles_C.append(counter_C)
        particles_D.append(counter_D)
        #velocities2_dt_total.append(velocities2_dt)
        #escrever_txt(velocities2_dt_total)
        time_pass.append(instant)

        #plt.plot(counter_A, counter_C, 'o', label=f"{instant}")

        screen.fill(background_color)
        
        if not paused_log:
            #Desenha todas as partículas na tela

            for particle in particles:
                particle.draw(screen)
            
                #print(velocities)
                
                #freq, bins = np.histogram(velocities, bins = 10, range = None)
                #freq = freq/len(particles)
                #print(freq)
                    
                #fig, axes = plt.subplots(1, 2)
                #plot_particles(axes[0], self.p_props, self.box_width, self.box_height, file_counter)
                #plot_v_histogram(axes[1], freq, bins, hist_x_limit = 600, hist_y_limit = 0.1)
                #print(plot_v_histogram)

                #f = [ self.maxwell_boltzmann(v, mass = self.mass, T = self.temperature) for v in vels]
                #axes[1].plot(vels, f, color='red')
                
                #plt.show()

                #fig.savefig(r'p_{}.png'.format(correct_counter(file_counter)) )
                #plt.close()
                    
                file_counter+=1


            pygame.display.flip()

            #particles_A.append(counter_A)
            #particles_B.append(counter_B)
            #particles_C.append(counter_C)
            #particles_D.append(counter_D)

            #plt.hist(velocities2_dt, bins=10)
            #plt.title(f"time {time}")
            #plt.show()
            
    #plt.plot(particles_A, particles_C)
    #plt.plot(counter_A, counter_C, 'o', label=f"{instant}")
     
    
    plt.plot(time_pass, particles_A, color="m") # (certo!!)
    plt.plot(time_pass, particles_C, color="b") # (certo!!)
    #plt.show() # (certo!!)

    modelo_fractionalA = Model(fractional, prefix="fractionalA_")
    print(f'parameter names: {modelo_fractionalA.param_names}')
    print(f'independent variables: {modelo_fractionalA.independent_vars}')
    params_fractionalA = modelo_fractionalA.make_params()
    
    modelo_squaredA = Model(squared, prefix="squaredA_")
    print(f'parameter names: {modelo_squaredA.param_names}')
    print(f'independent variables: {modelo_squaredA.independent_vars}')
    params_squaredA = modelo_squaredA.make_params()
        
    #modelo_squared2A = Model(squared_2, prefix="squared2A_")
    #print(f'parameter names: {modelo_squared2A.param_names}')
    #print(f'independent variables: {modelo_squared2A.independent_vars}')
    #params_squared2A = modelo_squared2A.make_params()
    
    modelo_firstorderA = Model(first_order, prefix="firstorderA_")
    print(f'parameter names: {modelo_firstorderA.param_names}')
    print(f'independent variables: {modelo_firstorderA.independent_vars}')
    params_firstorderA = modelo_firstorderA.make_params()

    params_fractionalA["fractionalA_C"].set(value=particles_A[0]+1, min=particles_A[0]-10, max=201)
    params_fractionalA["fractionalA_k"].set(value=-0.05, min=-1, max=5)
    
    params_squaredA["squaredA_C"].set(value=100, min=-100, max=201)
    params_squaredA["squaredA_k"].set(value=-10, min=-50, max=50)
    
    #params_squared2A["squared2A_CA"].set(value=particles_A[0]+1, min=particles_A[0]-1, max=100)
    #params_squared2A["squared2A_CB"].set(value=particles_B[0]+1, min=particles_B[0]-1, max=250)
    #params_squared2A["squared2A_B"].set(value=8, min=-100, max=100)
    #params_squared2A["squared2A_k"].set(value=-0.05, min=-50, max=50)
    
    params_firstorderA["firstorderA_C"].set(value=particles_A[0]+1, min=particles_A[0]-1, max=201)
    params_firstorderA["firstorderA_k"].set(value=-0.05, min=-1, max=5)
    
    #params_exponentialA["exponentialA_C"].set(value=particles_A[0]+1, min=particles_A[0]-1, max=particles_A[0]+5)
    #params_exponentialA["exponentialA_k"].set(value=-0.05, min=-1, max=5)

    resultado_fit_fractionalA = modelo_fractionalA.fit(particles_A, params_fractionalA, x=time_pass)
    print(resultado_fit_fractionalA.fit_report())
    
    print()
    
    resultado_fit_squaredA = modelo_squaredA.fit(particles_A, params_squaredA, x=time_pass)
    print(resultado_fit_squaredA.fit_report())
    
    print()
    
    #resultado_fit_squared2A = modelo_squared2A.fit(particles_A, params_squared2A, x=time_pass)
    #print(resultado_fit_squared2A.fit_report())
    
    print()
    
    resultado_fit_firstorderA = modelo_firstorderA.fit(particles_A, params_firstorderA, x=time_pass)
    print(resultado_fit_firstorderA.fit_report())
    
    resultado_fit_fractionalA.plot()
    resultado_fit_squaredA.plot()
    #resultado_fit_squared2A.plot()
    resultado_fit_firstorderA.plot()
    
    figure, axis = plt.subplots(2, 2, squeeze=False)
    
    fractional_fit = []
    squared_fit = []
    squared2_fit = []
    firstorder_fit = []
    for x in time_pass:
        fractional_fit.append(fractional(x, resultado_fit_fractionalA.best_values.get("fractionalA_C"), resultado_fit_fractionalA.best_values.get("fractionalA_k")))
        squared_fit.append(squared(x, resultado_fit_squaredA.best_values.get("squaredA_C"), resultado_fit_squaredA.best_values.get("squaredA_k")))
        #squared2_fit.append(squared_2(x, resultado_fit_squared2A.best_values.get("squared2A_CA"), resultado_fit_squared2A.best_values.get("squared2A_CB"), resultado_fit_squared2A.best_values.get("squared2A_B"), resultado_fit_squared2A.best_values.get("squared2A_k")))
        firstorder_fit.append(first_order(x, resultado_fit_firstorderA.best_values.get("firstorderA_C"), resultado_fit_firstorderA.best_values.get("firstorderA_k")))
    avg_x_vals_frac, deriv_vals_frac = derivada_xy(time_pass, fractional_fit)
    axis[0,0].plot(avg_x_vals_frac, deriv_vals_frac)
    axis[0,0].set_title("Derivada [A] - fit fractional")
    avg_x_vals_squ, deriv_vals_squ = derivada_xy(time_pass, squared_fit)
    axis[0,1].plot(avg_x_vals_squ, deriv_vals_squ)
    axis[0,1].set_title("Derivada [A] - fit squared")
    avg_x_vals_fo, deriv_vals_fo = derivada_xy(time_pass, firstorder_fit)
    axis[1,0].plot(avg_x_vals_fo, deriv_vals_fo)
    axis[1,0].set_title("Derivada [A] - fit first order")
    #avg_x_vals_squ2, deriv_vals_squ2 = derivada_xy(time_pass, squared2_fit)
    #axis[1,1].plot(avg_x_vals_squ2, deriv_vals_squ2)
    #axis[1,1].set_title("Derivada [A] - fit squared 2")
    plt.tight_layout()
    

    plt.show()
    

    pygame.quit()

    #print(particles_A)
    #print(particles_C)
    counter = 0
    #for veldt in velocities2_dt_total:
     #   plt.hist(veldt, bins=10, label=f"{counter}")
      #  plt.legend()
       # counter +=1
    #plt.show()


    
if __name__ == '__main__':
    main(options)
    
