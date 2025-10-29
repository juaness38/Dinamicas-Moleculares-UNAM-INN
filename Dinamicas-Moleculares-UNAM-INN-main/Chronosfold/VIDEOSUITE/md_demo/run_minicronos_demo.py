#!/usr/bin/env python3
"""Wrapper for MD + MiniChronoGPT + optional Milvus animation.

Delegates to scripts/minicronos_demo_video.py with PYTHONPATH set to chronosfold_scaffold/src.
"""
import argparse, os, subprocess, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--steps', type=int, default=400)
    ap.add_argument('--timestep', type=float, default=0.002)
    ap.add_argument('--temperature', type=float, default=300.0)
    ap.add_argument('--milvus', action='store_true')
    ap.add_argument('--collection', type=str, default='protein_chains')
    ap.add_argument('--limit', type=int, default=20)
    ap.add_argument('--pred-steps', type=int, default=60)
    ap.add_argument('--gif', type=str, default='reports/minicronos_demo.gif')
    ap.add_argument('--milvus-uri', type=str, default=None)
    ap.add_argument('--milvus-token', type=str, default=None)
    ap.add_argument('--fallback-venv', type=str, default=None, help='Path to python.exe with torch (e.g., .venv/Scripts/python.exe)')
    args = ap.parse_args()

    env = dict(os.environ)
    # Ensure chronosfold modules are importable
    env['PYTHONPATH'] = f"chronosfold_scaffold/src{os.pathsep}" + env.get('PYTHONPATH','')

    cmd = [sys.executable, 'scripts/minicronos_demo_video.py',
           '--steps', str(args.steps), '--timestep', str(args.timestep), '--temperature', str(args.temperature),
           '--pred-steps', str(args.pred_steps), '--gif', args.gif]
    if args.milvus:
        cmd.extend(['--milvus','--collection', args.collection, '--limit', str(args.limit)])
        if args.milvus_uri:
            cmd.extend(['--milvus-uri', args.milvus_uri])
        if args.milvus_token:
            cmd.extend(['--milvus-token', args.milvus_token])
    if args.fallback_venv:
        cmd.extend(['--fallback-venv', args.fallback_venv])

    os.makedirs(os.path.dirname(args.gif) or '.', exist_ok=True)
    print('Running:', ' '.join(cmd))
    subprocess.check_call(cmd, env=env)

if __name__ == '__main__':
    main()
