#!/usr/bin/env python3
"""Wrapper to generate a simple 2D synthetic MP4 using ChronosFold plotting utilities.

Delegates to chronosfold_scaffold FastAPI router function animate_mp4().
"""
import argparse, asyncio, base64, pathlib, sys
sys.path.insert(0, 'chronosfold_scaffold/src')
from chronosfold.api.router import animate_mp4 as _animate_mp4

async def _run(out: str, frames: int, fps: int, color_by: str, codec: str):
    resp = await _animate_mp4(frames=frames, fps=fps, color_by=color_by, codec=codec)
    data = base64.b64decode(resp['mp4_base64'])
    p = pathlib.Path(out)
    p.write_bytes(data)
    print(f"Wrote {p} frames={resp['frames']} fps={resp['fps']} codec={resp['codec']}")

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=str, default='alanine_dipeptide_demo.mp4')
    ap.add_argument('--frames', type=int, default=200)
    ap.add_argument('--fps', type=int, default=20)
    ap.add_argument('--color-by', type=str, default='density', choices=['density','cluster'])
    ap.add_argument('--codec', type=str, default='libx264')
    args = ap.parse_args()
    asyncio.run(_run(args.out, args.frames, args.fps, args.color_by, args.codec))
