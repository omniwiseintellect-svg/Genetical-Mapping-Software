"""Generate a small artificial FASTA example using seed 20261027."""
from pathlib import Path
import random

SEED = 20261027

def generate():
    rng = random.Random(SEED)
    base = ''.join(rng.choices('ACGT', k=120))
    lines = []
    for population, offset in [('DemoA', 0), ('DemoB', 3), ('DemoC', 21), ('DemoD', 47)]:
        for index in range(2):
            sequence = list(base)
            for position in range(offset, offset + 12):
                site = (position + index * 17) % len(sequence)
                sequence[site] = rng.choice([base for base in 'ACGT' if base != sequence[site]])
            lines.extend([
                f'>{population}_sequence{index + 1} type:complete source:synthetic seed:{SEED}',
                ''.join(sequence),
            ])
    return '\n'.join(lines) + '\n'

if __name__ == '__main__':
    target = Path(__file__).with_name('examples') / 'synthetic_demo.fasta'
    target.parent.mkdir(exist_ok=True)
    target.write_text(generate(), encoding='utf-8', newline='\n')
    print(target.name)
