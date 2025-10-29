#!/usr/bin/env python3
"""Minimal client to call /chronosfold/animate_mp4_raw and save the video.
"""
import argparse, base64, json, sys
import urllib.request

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--host', type=str, default='http://localhost:8000')
    ap.add_argument('--out', type=str, default='demo.mp4')
    ap.add_argument('--frames', type=int, default=240)
    ap.add_argument('--fps', type=int, default=24)
    ap.add_argument('--color-by', type=str, default='density')
    ap.add_argument('--codec', type=str, default='libx264')
    args = ap.parse_args()

    import urllib.parse
    q = urllib.parse.urlencode({'frames': args.frames, 'fps': args.fps, 'color_by': args.color_by, 'codec': args.codec})
    url = f"{args.host}/chronosfold/animate_mp4_raw?{q}"

    print('GET', url)
    with urllib.request.urlopen(url) as resp:
        data = resp.read()
    with open(args.out, 'wb') as f:
        f.write(data)
    print('Saved to', args.out)

if __name__ == '__main__':
    main()
