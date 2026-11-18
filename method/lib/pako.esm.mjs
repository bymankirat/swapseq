// Small stand-in for pako used by the bundled jsfive: inflate/deflate via Node's built-in zlib.
import zlib from 'zlib';
export function inflate(buf){return new Uint8Array(zlib.inflateSync(Buffer.from(buf.buffer?new Uint8Array(buf.buffer,buf.byteOffset,buf.byteLength):buf)))}
export function deflate(buf){return new Uint8Array(zlib.deflateSync(Buffer.from(buf)))}
