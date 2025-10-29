#!/usr/bin/env python3
"""Wrapper for PyMOL-based semantic trajectory video.

Runs: python -m ste.render.pymol_video with given args.
"""
import argparse, os, subprocess, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pdb', required=True)
    ap.add_argument('--traj', required=False, default=None)
    ap.add_argument('--annotations', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--fps', type=int, default=30)
    ap.add_argument('--color-mode', choices=['state','novelty'], default='state')
    ap.add_argument('--width', type=int, default=800)
    ap.add_argument('--height', type=int, default=600)
    args = ap.parse_args()

    env = dict(os.environ)
    env['PYTHONPATH'] = f"{env.get('PYTHONPATH','')}{os.pathsep}chronosfold_scaffold/src"

    cmd = [sys.executable, '-m', 'ste.render.pymol_video',
           '--pdb', args.pdb, '--annotations', args.annotations,
           '--out', args.out, '--fps', str(args.fps), '--color-mode', args.color_mode,
           '--width', str(args.width), '--height', str(args.height)]
    if args.traj:
        cmd.extend(['--traj', args.traj])

    os.makedirs(args.out, exist_ok=True)
    print('Running:', ' '.join(cmd))
    subprocess.check_call(cmd, env=env)

if __name__ == '__main__':
    main()
