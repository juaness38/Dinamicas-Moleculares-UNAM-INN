import base64, pathlib
from chronosfold.data.loaders import generate_synthetic_frames
from chronosfold.pipeline.semanticize import build_embeddings
from chronosfold.visualization.plotting import trajectory_animation_mp4

path = pathlib.Path('videos/VIDEOSUITE/2d/manual_fallback.mp4')
path.parent.mkdir(parents=True, exist_ok=True)
frames = generate_synthetic_frames(n_frames=80, seed=12)
traj = build_embeddings(frames)
b64 = trajectory_animation_mp4(traj.to_json(), color_by='density', fps=12, codec='libx264')
path.write_bytes(base64.b64decode(b64))
print('wrote', path)
