#!/usr/bin/env python3
"""Wrapper to generate a simple 2D synthetic MP4 using ChronosFold plotting utilities.

Default path uses plotting helpers directly (fast, no heavy imports). Optionally can
delegate to FastAPI router with --use-router.
"""
import argparse, asyncio, base64, pathlib, sys, os
sys.path.insert(0, 'chronosfold_scaffold/src')
print('[2d] wrapper starting', flush=True)

async def _run(out: str, frames: int, fps: int, color_by: str, codec: str, use_router: bool = False):
    p = pathlib.Path(out)
    p.parent.mkdir(parents=True, exist_ok=True)
    if use_router:
        print('[2d] using router path', flush=True)
        from chronosfold.api.router import animate_mp4 as _animate_mp4  # type: ignore
        resp = await _animate_mp4(frames=frames, fps=fps, color_by=color_by, codec=codec)
        data = base64.b64decode(resp['mp4_base64'])
        p.write_bytes(data)
        print(f"Wrote {p} frames={resp['frames']} fps={resp['fps']} codec={resp['codec']}")
    else:
        print('[2d] using direct plotting path', flush=True)
        from chronosfold.pipeline.semanticize import build_embeddings  # type: ignore
        from chronosfold.data.loaders import generate_synthetic_frames  # type: ignore
        from chronosfold.visualization.plotting import trajectory_animation_mp4  # type: ignore
        f = generate_synthetic_frames(n_frames=frames, seed=12)
        traj = build_embeddings(f)
        b64 = trajectory_animation_mp4(traj.to_json(), color_by=color_by, fps=fps, codec=codec)
        p.write_bytes(base64.b64decode(b64))
        print(f"Wrote {p} frames={frames} fps={fps} codec={codec} (fallback)")

if __name__ == '__main__':
    try:
        (pathlib.Path('videos/VIDEOSUITE/2d')).mkdir(parents=True, exist_ok=True)
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=str, default='alanine_dipeptide_demo.mp4')
    ap.add_argument('--frames', type=int, default=200)
    ap.add_argument('--fps', type=int, default=20)
    ap.add_argument('--color-by', type=str, default='density', choices=['density','cluster'])
    ap.add_argument('--codec', type=str, default='libx264')
    ap.add_argument('--use-router', action='store_true', help='Use FastAPI router animate_mp4 (optional)')
    args = ap.parse_args()
    try:
        print(f"[2d] invoking _run out={args.out} frames={args.frames} fps={args.fps} color_by={args.color_by} codec={args.codec} use_router={args.use_router}", flush=True)
        asyncio.run(_run(args.out, args.frames, args.fps, args.color_by, args.codec, use_router=args.use_router))
    except Exception as e:
        # Write a minimal debug log
        dbg = pathlib.Path('videos/VIDEOSUITE/2d/runner_debug.log')
        try:
            dbg.write_text(f"error: {e}\n", encoding='utf-8')
        except Exception:
            pass
        raise
