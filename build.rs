// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=assets/itme.ico");

    if env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows") {
        return;
    }

    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR"));
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR"));
    let ico_path = manifest_dir.join("assets").join("itme.ico");
    if !ico_path.exists() {
        println!(
            "cargo:warning=Missing icon file at {}, skipping Windows resource embedding",
            ico_path.display()
        );
        return;
    }

    let ico_escaped = ico_path.to_string_lossy().replace('\\', "/");
    let pkg_version = env::var("CARGO_PKG_VERSION").unwrap_or_else(|_| "2.3.0".to_string());
    let parts: Vec<u16> = pkg_version
        .split('.')
        .filter_map(|p| p.split('-').next()?.parse::<u16>().ok())
        .collect();
    let v0 = *parts.first().unwrap_or(&2);
    let v1 = *parts.get(1).unwrap_or(&3);
    let v2 = *parts.get(2).unwrap_or(&0);
    let v3 = *parts.get(3).unwrap_or(&0);

    let rc_content = format!(
        r#"1 ICON "{ico_escaped}"

1 VERSIONINFO
FILEVERSION {v0},{v1},{v2},{v3}
PRODUCTVERSION {v0},{v1},{v2},{v3}
FILEFLAGSMASK 0x3fL
FILEFLAGS 0x0L
FILEOS 0x00040004L
FILETYPE 0x00000001L
FILESUBTYPE 0x0L
BEGIN
    BLOCK "StringFileInfo"
    BEGIN
        BLOCK "040904B0"
        BEGIN
            VALUE "CompanyName", "CRV Science Network\0"
            VALUE "FileDescription", "ITME — Information-Theoretic Model of Evolution\0"
            VALUE "FileVersion", "{pkg_version}\0"
            VALUE "InternalName", "itme_evolution_sim\0"
            VALUE "LegalCopyright", "Copyright CRV Science Network (GPL-3.0-or-later)\0"
            VALUE "OriginalFilename", "4_ITME.exe\0"
            VALUE "ProductName", "ITME\0"
            VALUE "ProductVersion", "{pkg_version}\0"
        END
    END
    BLOCK "VarFileInfo"
    BEGIN
        VALUE "Translation", 0x0409, 1200
    END
END
"#
    );

    let rc_path = out_dir.join("itme.rc");
    if let Err(e) = fs::write(&rc_path, rc_content) {
        println!("cargo:warning=Failed to write {}: {e}", rc_path.display());
        return;
    }

    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();
    if target_env == "msvc" {
        let res_path = out_dir.join("app_res.res");
        if let Some(rc_exe) = find_rc_exe() {
            let status = Command::new(&rc_exe)
                .args(["/nologo", "/fo"])
                .arg(&res_path)
                .arg(&rc_path)
                .status();
            match status {
                Ok(s) if s.success() && res_path.exists() => {
                    println!("cargo:rustc-link-arg={}", res_path.display());
                }
                Ok(s) => {
                    println!("cargo:warning=rc.exe exited with status {s}, skipping resource embed");
                }
                Err(e) => {
                    println!("cargo:warning=Failed to run rc.exe ({e}), skipping resource embed");
                }
            }
        } else {
            println!("cargo:warning=rc.exe not found for MSVC target, skipping resource embed");
        }
    } else {
        let obj_path = out_dir.join("app_res.o");
        let windres_candidates = [
            "windres",
            "x86_64-w64-mingw32-windres",
            "windres.exe",
            "x86_64-w64-mingw32-windres.exe",
        ];
        let mut compiled = false;
        for tool in windres_candidates {
            if let Ok(status) = Command::new(tool)
                .arg("-i")
                .arg(&rc_path)
                .arg("-O")
                .arg("coff")
                .arg("-o")
                .arg(&obj_path)
                .status()
            {
                if status.success() && obj_path.exists() {
                    println!("cargo:rustc-link-arg={}", obj_path.display());
                    compiled = true;
                    break;
                }
            }
        }
        if !compiled {
            println!("cargo:warning=windres not found or failed, skipping Windows resource embed");
        }
    }
}

fn find_rc_exe() -> Option<PathBuf> {
    if let Ok(status) = Command::new("rc.exe").arg("/?").output() {
        if !status.stdout.is_empty() || !status.stderr.is_empty() {
            return Some(PathBuf::from("rc.exe"));
        }
    }
    let kits_root = Path::new(r"C:\Program Files (x86)\Windows Kits\10\bin");
    if kits_root.exists() {
        if let Ok(entries) = fs::read_dir(kits_root) {
            let mut versions: Vec<PathBuf> = entries
                .filter_map(|e| e.ok().map(|x| x.path()))
                .filter(|p| p.is_dir())
                .collect();
            versions.sort();
            versions.reverse();
            for ver in versions {
                let candidate = ver.join("x64").join("rc.exe");
                if candidate.exists() {
                    return Some(candidate);
                }
            }
        }
    }
    None
}
