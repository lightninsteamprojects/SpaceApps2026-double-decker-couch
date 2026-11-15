-- Database setup for SPHEREx Sky Blink (Supabase / Postgres).
-- Run this once in your own project: Supabase dashboard -> SQL Editor -> paste -> Run.
-- It is safe to re-run on an empty project; on an existing one, read it first (it creates tables and policies).
--
-- No keys or secrets are in this file. The page talks to Supabase with the project's public ("publishable" / anon)
-- key, so Row Level Security (RLS) below is what actually limits what visitors can do.


-- ===========================================================================================================
-- comments: the Discussion tab. Anyone can post; a comment only becomes visible once YOU set approved = true
-- (Table Editor -> comments -> tick "approved"). That is the moderation step.
-- ===========================================================================================================
create table if not exists public.comments (
  id         bigint generated always as identity primary key,
  created_at timestamptz not null default now(),
  field_id   text not null check (char_length(field_id) between 1 and 60),
  band       text check (char_length(band) <= 10),
  nickname   text check (char_length(nickname) <= 30),
  body       text not null check (char_length(body) between 1 and 500),
  approved   boolean not null default false
);

alter table public.comments enable row level security;

grant insert, select on public.comments to anon;

drop policy if exists "anyone can post" on public.comments;
create policy "anyone can post" on public.comments
  for insert to anon
  with check (approved = false);            -- visitors cannot post pre-approved comments

drop policy if exists "anyone can read approved" on public.comments;
create policy "anyone can read approved" on public.comments
  for select to anon
  using (approved = true);                  -- unapproved comments stay hidden from the public


-- ===========================================================================================================
-- classifications: the answers behind "What did you see?" (one row per answer).
-- CHECK AGAINST THE LIVE TABLE BEFORE USE: this was derived from what index.html inserts
-- (target_id, answer, session_id, comment, band), not read from the team's real database, so the real table
-- may have extra columns, different types or constraints.
-- ===========================================================================================================
create table if not exists public.classifications (
  id         bigint generated always as identity primary key,
  created_at timestamptz not null default now(),
  target_id  text not null,      -- field id, e.g. 'pleiades' (the id in targets.json)
  answer     text not null,      -- 'moved' | 'brightness' | 'nothing'
  session_id text,               -- random id the browser keeps in localStorage ('sid'); not a login
  comment    text,               -- optional note (the page trims it to 500 characters)
  band       text                -- e.g. 'SPHEREx-D2'; null for fields without band data
);

alter table public.classifications enable row level security;

grant insert on public.classifications to anon;          -- insert only: visitors can't read other people's answers

drop policy if exists "anyone can submit an answer" on public.classifications;
create policy "anyone can submit an answer" on public.classifications
  for insert to anon
  with check (true);
