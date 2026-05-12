#!/usr/bin/env node
/**
 * Headless screenshot utility for Moonbase.
 *
 * Usage:
 *   node scripts/screenshot.mjs [options]
 *
 * Options:
 *   --url        Page URL             (default: http://localhost:5177)
 *   --out        Output PNG path      (default: /tmp/moonbase-screenshot.png)
 *   --width      Viewport width       (default: 1920)
 *   --height     Viewport height      (default: 1080)
 *   --wait       Wait ms after load   (default: 4000)
 *   --click      CSS selectors to click, comma-separated (executed in order)
 *   --clickwait  Wait ms after each click (default: 2000)
 *   --base       Switch to base before screenshot (moon/mars/titan/etc, default: moon)
 *   --full       Full-page screenshot (default: false)
 *
 * Exit codes:
 *   0  Screenshot taken, no errors
 *   1  Fatal error (couldn't load page, browser crash, etc.)
 *   2  Screenshot taken but page had errors (R3F, runtime exceptions, Vite overlay)
 *
 * Examples:
 *   node scripts/screenshot.mjs
 *   node scripts/screenshot.mjs --base titan
 *   node scripts/screenshot.mjs --width 390 --height 844
 */

import { chromium } from 'playwright'

// Errors matching these patterns are harmless HMR/Zustand noise — ignore them
const IGNORED_PATTERNS = [
  /structuredClone/,
  /Zustand/i,
]

function isIgnored(msg) {
  return IGNORED_PATTERNS.some(p => p.test(msg))
}

// ─── Parse args ──────────────────────────────────────────────────────────────

function parseArgs() {
  const args = process.argv.slice(2)
  const opts = {
    url: 'http://localhost:5177',
    out: '/tmp/moonbase-screenshot.png',
    width: 1920,
    height: 1080,
    wait: 4000,
    click: '',
    clickwait: 2000,
    base: 'moon',
    full: false,
  }

  for (let i = 0; i < args.length; i++) {
    const key = args[i].replace(/^--/, '')
    if (key === 'full') {
      opts.full = true
      continue
    }
    const val = args[++i]
    if (key in opts) {
      opts[key] = ['width', 'height', 'wait', 'clickwait'].includes(key)
        ? parseInt(val, 10)
        : val
    }
  }
  return opts
}

// ─── Main ────────────────────────────────────────────────────────────────────

async function main() {
  const opts = parseArgs()
  const errors = []

  const browser = await chromium.launch({
    headless: true,
    args: [
      '--enable-webgl',
      '--enable-webgl2',
      '--use-gl=angle',
      '--use-angle=swiftshader',
      '--enable-features=Vulkan',
    ],
  })
  const context = await browser.newContext({
    viewport: { width: opts.width, height: opts.height },
    deviceScaleFactor: 2, // retina-quality screenshots
  })
  const page = await context.newPage()

  // Collect errors
  page.on('pageerror', e => {
    if (!isIgnored(e.message)) {
      errors.push(`PAGE ERROR: ${e.message}`)
      console.error(`PAGE ERROR: ${e.message}`)
    }
  })
  page.on('console', msg => {
    if (msg.type() === 'error') {
      const text = msg.text()
      if (!isIgnored(text)) {
        errors.push(`CONSOLE ERROR: ${text}`)
        console.error(`CONSOLE ERROR: ${text}`)
      }
    }
  })

  try {
    await page.goto(opts.url, { waitUntil: 'networkidle', timeout: 15000 })
  } catch {
    // networkidle can be flaky with animation loops — fall back to load
    await page.goto(opts.url, { waitUntil: 'load', timeout: 15000 })
  }

  // Check for Vite error overlay (import failures, syntax errors, etc.)
  const hasViteOverlay = await page.$('vite-error-overlay').then(el => !!el).catch(() => false)
  if (hasViteOverlay) {
    const overlayText = await page.evaluate(() => {
      const el = document.querySelector('vite-error-overlay')
      return el?.shadowRoot?.textContent?.slice(0, 500) || 'Vite error overlay detected'
    }).catch(() => 'Vite error overlay detected')
    errors.push(`VITE OVERLAY: ${overlayText}`)
    console.error(`VITE OVERLAY: ${overlayText}`)
  }

  // Switch base if requested via hover dropdown
  if (opts.base && opts.base !== 'moon') {
    // Hover the base switcher to open dropdown
    await page.hover('.group .font-mono.text-cyan')
    await page.waitForTimeout(300)
    // Click the matching dropdown button — find by partial text in button list
    const buttons = await page.$$('.group.relative button')
    let switched = false
    for (const btn of buttons) {
      const text = await btn.textContent()
      const nameMap = { mars: 'ARES', titan: 'HUYGENS', europa: 'EUROPA', ceres: 'DAWN', venus: 'APHRODITE', phobos: 'PHOBOS' }
      const keyword = nameMap[opts.base] || opts.base.toUpperCase()
      if (text && text.includes(keyword)) {
        await btn.click()
        await page.waitForTimeout(opts.wait)
        switched = true
        break
      }
    }
    if (!switched) {
      errors.push(`BASE SWITCH: Could not find button for base "${opts.base}"`)
      console.error(`BASE SWITCH: Could not find button for base "${opts.base}"`)
    }

    // Check for Vite overlay again after base switch (new code paths loaded)
    const hasOverlayAfter = await page.$('vite-error-overlay').then(el => !!el).catch(() => false)
    if (hasOverlayAfter && !hasViteOverlay) {
      const overlayText = await page.evaluate(() => {
        const el = document.querySelector('vite-error-overlay')
        return el?.shadowRoot?.textContent?.slice(0, 500) || 'Vite error overlay detected'
      }).catch(() => 'Vite error overlay after base switch')
      errors.push(`VITE OVERLAY (after switch): ${overlayText}`)
      console.error(`VITE OVERLAY (after switch): ${overlayText}`)
    }
  }

  // Wait for Three.js to render
  await page.waitForTimeout(opts.wait)

  // Execute click sequence if provided
  if (opts.click) {
    const selectors = opts.click.split(',').map(s => s.trim()).filter(Boolean)
    for (const sel of selectors) {
      try {
        await page.click(sel, { timeout: 5000, force: true })
        await page.waitForTimeout(opts.clickwait)
      } catch (e) {
        console.error(`Warning: click "${sel}" failed: ${e.message}`)
      }
    }
  }

  // Take screenshot
  await page.screenshot({
    path: opts.out,
    fullPage: opts.full,
    timeout: 60000,
  })

  await browser.close()

  // Report results
  console.log(opts.out)

  if (errors.length > 0) {
    console.error(`\n⚠ ${errors.length} error(s) detected:`)
    for (const err of errors) {
      console.error(`  • ${err.slice(0, 200)}`)
    }
    process.exit(2)
  }
}

main().catch(e => {
  console.error(e.message)
  process.exit(1)
})
