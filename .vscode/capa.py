import pygame
import sys
import math

# Inicializa o Pygame
pygame.init()

# Configurações da tela
width, height = 1584, 396
screen = pygame.display.set_mode((width, height))
pygame.display.set_caption("LinkedIn Cover")

# Cores e fonte
WHITE = (255, 255, 255)
font = pygame.font.SysFont('arial', 64, bold=True)

# Relógio para animação
clock = pygame.time.Clock()

running = True
t = 0

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    # Gradiente animado de fundo
    t += 0.02
    for y in range(height):
        r = int(100 + 50 * math.sin(t))
        g = int(100 + 50 * math.cos(t + 2))
        b = int(100 + 50 * math.sin(t + 4))
        alpha = int(255 - (y / height) * 105)
        pygame.draw.line(screen, (r, g, b), (0, y), (width, y))

    scale_factor = 1 + 0.05 * math.sin(t * 0.2)
    text = font.render("Julia Cristina", True, (*WHITE, 200))
    text_rect = text.get_rect(center=(width / 2, height / 2))

    scaled_text = pygame.transform.smoothscale(text, 
        (int(text_rect.width * scale_factor), int(text_rect.height * scale_factor)))
    scaled_rect = scaled_text.get_rect(center=(width / 2, height / 2))
    
    screen.blit(scaled_text, scaled_rect)

    pygame.display.flip()
    clock.tick(60)

    # Salvar a imagem (pressione 'S' para salvar)
    if pygame.key.get_pressed()[pygame.K_s]:
        pygame.image.save(screen, "linkedin_cover.png")
        print("Imagem salva como linkedin_cover.png")

pygame.quit()
sys.exit()