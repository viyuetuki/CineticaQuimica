class Event(object):
    
    """Uma classe para o evento de colisão que contém informações sobre o tipo, tempo e partículas envolvidas."""
    
    def __init__(self,ct,t,p1_index,p2_index=None,wall=None):

        """Colisões, que podem ser tanto entre partículas ou com as paredes.
        
        Args:
            ct: tipo da colisão ("w" é uma colisão com a parede, qualquer outra coisa é entre duas partículas)
            t: tempo da colisão
            p1_index: parâmetro de localização de uma primeira partícula
            p2_index: parâmetro de localização de uma segunda partícula
            (opcional)
            wall: nome da parede envolvida na colisão
        """
        
        # Atributos do evento
        self.wc_log = False # is it wall collision? (otherwise binary particle col)
        self.t = t
        self.p1_index = p1_index
        self.p2_index = p2_index
        self.wall = wall   
        
        #Configura os atributos do evento de acordo com o tipo da colisão
        if ct == "w":
            self.wc_log = True
        else:
            # Adverte caso o index da segunda partícula não foi passado
            if p2_index == None:
                raise RuntimeWarning("Warning: Event: second particle index undefined")